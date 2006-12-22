#!/usr/local/bin/perl
#-T
# use strict;
$|++;

##### program to calculate length of terminal restriction fragments 
##### for many enzymes and compare results among groups
##### input is a multiple sequence alignment in fasta format,
##### groups are coded as the first two characters of the taxa name
##### and tab delimited file with enzyme and cut site on each line
##### get rid of whitespace in sequences names and start sequences at beginning of primer
##### to run from the command line: perl repk_[version].pl [alignment file] [enzyme file]
##### to run by double-clicking: rename alignment file to "alignment.fasta" and enzyme file to enzymes.txt and run.

##### in _1, put reverse complementation inside loop, added variable
##### in _2, fixed bug in stringency setting
##### in _2z1, prints all fragment lengths now
##### todo in _3, add ability to use custom cut sites (e.g. uncut)
##### todo in _i, add ability to set stringency for number of allowable misses per group

#use strict;
# use CGI::Carp qw(fatalsToBrowser);
use File::Copy;	#a module necessary for copying files

######################
######################
#####################


# First, get the CGI variables into a list of strings
%cgivars= &getcgivars ;

# Print the CGI response header, required for all HTML output
# Note the extra \n, to send the blank line
print "Content-type: text/html\n\n" ;

# Finally, print out the complete HTML response page
print <<EOF ;
<html>
<head><title>Repicker</title></head>
<body>
<h2>Results Page</h2>
The script is not finished running until it says <h3>STOP</h3>
<br>Open links into new tabs or save them to disk, but don't click directly on them or you'll cancel the run.
<br>It might take a while, especially if it has been used frequently and given a lower priority by the server.
<br>There is a downloadable script you can run locally on your own computer, go back to the main page to get it.
<br><br>
EOF

# Print the CGI variables sent by the user.
# Note that the order of variables is unpredictable.
# Also note this simple example assumes all input fields had unique names,
#   though the &getcgivars() routine correctly handles similarly named
#   fields-- it delimits the multiple values with the \0 character, within 
#   $cgivars{$_}.
# foreach (keys %cgivars) {
#    print "<li>[$_] = [$cgivars{$_}]\n" ;
# }


##################
###################
#################


my($sec,$min,$hour,$day,$mon,$year,undef,undef,undef) = localtime time;$year += 1900;
my $timestring = "REPK_"."$year$mon$day-$hour$min-$sec"; 
mkdir $timestring, 0755;

##### variables that can be set
my $cutoff = $cgivars{'cutoff'}; #5;	##### furthest apart two fragments can be in length and still be considered the same fragments
				##### could split up the cutoff into 2: one for INTERgroup, another for INTRAgroup
my $stringency = $cgivars{'stringency'}/100; #0.0; 	# an enzyme must distinguish MORE than this percentage of groups to be acceptable
#my $matrixwarning = $cgivars{'matrixwarning'}; #0.1;# doesn't work, but should: display warning about groups that are cut by fewer than this fraction of enzymes
my $bpcutoff_low = $cgivars{'bpcutoff_low'}; #75;	# shortest acceptable fragment
my $bpcutoff_high = $cgivars{'bpcutoff_high'}; #900;# longest acceptable fragment
my $matchlimit = $cgivars{'matchlimit'}; #100;	# the most number of matches we want printed
my $splitter = $cgivars{'splitter'}; #"_";
my $splitnum = $cgivars{'splitnum'} -1;


##### use default filenames if none are provided on the command line
# @ARGV[0] = "alignment.fasta" if length @ARGV[0] == 0;
@ARGV[1] = "enzymes_comm_iso.txt" if length @ARGV[1] == 0;

##### input filenames
# my $alignmentFile = @ARGV[0];			#"SequenceInputFile.fasta";
my $enzymeFile = @ARGV[1];			#"EnzymeInputFile.txt";

##### sets are (0,4),(5,5),(10,6),(15,6)
# my $taxabid = 6;		#where group names begins (0 is first)
# my $taxaeid = 6;		#how many characters is group name


##### this program makes copies of the enzyme file, sequence alignment file, and itself to the results directory, for later consultation
#copy("$enzymeFile","./"."$timestring"."/");
#copy("$alignmentFile","./"."$timestring"."/");
#copy("$0","./"."$timestring"."/");	#$0 is the variable that holds the perl script's name

# ##### output filenames
# ##if you comment out use strict you can here comment out any output files you don't want printed
my $yescutsfile = "./$timestring/yes_cuts.txt";		# enzymes that cut and meet stringency
my $successfile = "./$timestring/success.txt";	# enzymes that distinguish all the groups
my $failfile = "./$timestring/no_cuts.txt";		# enzymes that do not cut/distinguish any groups
my $matrixfile = "./$timestring/enzmatrix.txt";	# matrix of group combinations each passing enzyme cuts
my $outfile = "./$timestring/finalout.txt";		# list of 4-enzyme groups that will identify all groups
my $fragfile = "./$timestring/fragfile.txt";		# fragment sizes for each passing enzyme


##### my $fragprintall = 1;	# set to 1 if all fragments are to be printed, not just the best, NOT YET INSTITUTED

die "CUTOFF can not be greater than or equal to BPCUTOFF_LOW or trouble will ensue" if ($bpcutoff_low <= $cutoff);

##### widely used hashes
my %success;
my %fragments;
my %fragsave;
my %list;
my %partlist;
my %cutlength;
my %matchlength;
my %reverseEnzymes;
my %isoschiz;
my %goodcombo;

##### variables that are printed out of scope
my $listcount;
my $partcount;
my @nocutters;
my $combonumber=0;
my %finalsave;
my @subs;
# my @zerogroup;

#######################    
####### PART 1 ########  In which the fragments are cut with the enzymes, the terminal fragment lengths are saved,
#######################  the groups are differentiated, evil is vanquished, and peace is brought to the entire galaxy.

my %sequences = readFASTAfile (%cgivars); #cgi $alignmentFile);
@subs = checkGroups (%sequences);

my $numgroups = @subs;
print "numgroups = $numgroups\n<br>";
my $numcombos = $numgroups*($numgroups-1)/2;
print "numcombos = $numcombos\n<br>";
my $maxfail = $numcombos - int($numcombos*$stringency);
print "maxfail = $maxfail\n<br>";

my @grpgrp = groupGroups (@subs);
my @revgrpgrp = reverse @grpgrp;

my %enzymes = readEnzymeFile ($enzymeFile);
print((scalar keys %enzymes) . " enzymes read from <a href=\"$enzymeFile\">$enzymeFile</a>\n<br>");
#   Groups distinguished by fewer than " . (scalar keys %enzymes) * $matrixwarning . " enzymes will be printed to screen.\n<br>");

my $enzhashsize = keys(%enzymes);
ENZYME:	foreach my $RE (sort (keys %enzymes)) {
# 	print "testing $RE\n";  #####\t$enzymes{$RE}\t$reverseEnzymes{$RE}\t$cutlength{$RE}\n";
# 	print "testing $RE\tISO: @{$isoschiz{$RE}}\n";	##### debugging
# 	print (($enzhashsize--)."...");
	$enzhashsize--;
	if (($enzhashsize % 10) == 0) {print ("$enzhashsize")}
	else {print "."}

	##### reverse complements each sequence to search for opposite strand cutsites
	$reverseEnzymes{$RE} = reverse($enzymes{$RE});
	$reverseEnzymes{$RE} =~ tr/ACGT[]/TGCA][/;
	
	my $site = $enzymes{$RE};
	my $revsite = $reverseEnzymes{$RE};
	
	my %seen = ();
	my %duplicates = ();
	my $failcount = 0;


	##### finds fragment lengths for all seqs, sets to negative value if not within bpcutoffs
	foreach my $taxa (sort (keys %sequences)) {
		$fragments{$taxa} = ();
		my $terminalFragment = 0;
		my $tempfrag = 0;
		my $forwardFragment;
		my $reverseFragment;
		
		##### checks the sequence for the forward recognition sequence, saves it to terminalFragment if passes stringency, tempfrag if not
		if ($sequences{$taxa} =~ m/($site){1}/i) {
			$forwardFragment = length ($`) + $cutlength{$RE};
			if (($forwardFragment <$bpcutoff_high) and ($forwardFragment >$bpcutoff_low)) {
				$terminalFragment = $forwardFragment;
			}
			else {
				$tempfrag = $forwardFragment;
			}
			
		}
		##### checks the sequence for the reverse recognition sequence, saves it to terminalFragment if passes stringency 
		##### and is shorter than forward terminalFragment, tempfrag if it is shorter but doesn't pass stringency
		if ($sequences{$taxa} =~ m/($revsite){1}/i) {
			$reverseFragment = length ($`) + $cutlength{$RE};
			if (($reverseFragment <$bpcutoff_high) and ($reverseFragment >$bpcutoff_low) and ($reverseFragment < $forwardFragment)) {
				$terminalFragment = $reverseFragment;	
			}
			else {
				($tempfrag = $reverseFragment) if ($reverseFragment < $tempfrag);
			}
		}
		
		##### saves the fragment length that passes stringency for later printing and analysis
		if ($terminalFragment) {
			$fragments{$taxa} = "$terminalFragment";
			push @{$fragsave{$RE}}, "$terminalFragment";
	}
		##### saves the fragment length that does not pass stringency for printing later
		elsif ($tempfrag) {
			$fragments{$taxa} = 0;
			push @{$fragsave{$RE}},1*$tempfrag; #was -1
		}
		##### saves the entire fragment length because there were no cuts
		else {
			$fragments{$taxa} = 0;
			push @{$fragsave{$RE}},1*length ($sequences{$taxa}); #was -1
		}

	}
	
	##### saves any enzymes that cuts no sequences or cuts them all alike for later print out to failure file
	my %flip_frag = reverse %fragments;
	unless (keys %flip_frag > 1) {
		push @nocutters, $RE;
		next ENZYME;
		}

	##### compares groups, finds which groups have similar patterns
	CHECK:	foreach my $checkkey (sort (keys %fragments)) {  ##### first keys are sequence names
# 		my $set_1 = substr ($checkkey, $taxabid, $taxaeid);  ##### pulls group name out of first sequence name
		my @set_1array = split (/$splitter/, $checkkey);
		my $set_1 = @set_1array[$splitnum];
		
		foreach my $fragtaxa (sort (keys %fragments)) {	##### second keys are sequence names
			next if ($fragtaxa le $checkkey);			##### don't check the same taxa against itself or any less than it
# 			my $set_2 = substr ($fragtaxa, $taxabid, $taxaeid);	##### pulls group name out of second sequence name
			my @set_2array = split (/$splitter/, $fragtaxa);
			my $set_2 = @set_2array[$splitnum];
			next if ($seen{$set_1} =~ m/($set_2)/);		##### don't check again if it's already been flagged
			my $lengthDifference =  abs ($fragments{$checkkey} - $fragments{$fragtaxa});	##### find the difference in length between fragments
			
# 			print "$RE\t$set_1\t$set_2\t$lengthDifference\n";	##### debugging
			
			################# if you want to be more stringent so that all members of the same group must have the same length ######
# 			if (($lengthDifference >= $cutoff) and ($set_1 eq $set_2)) {
# 				push @nocutters, $RE;
# 				next ENZYME;}
			#########################################################################################################################
			
			next if ($set_1 eq $set_2);				# we don't care as long as they're not different
			next if ($lengthDifference > $cutoff);	# we DO want different groups to have different lengths
# 			print "$set_1 and $set_2 don't match: difference of $lengthDifference is less than cutoff of $cutoff\n";	##### debugging
			$seen{$set_1} .= "$set_2\t";			# flag the groups to know they won't match
# 			print "$set_1: $seen{$set_1}\n";		##### debugging
		}
	}

	##### saves which sets failed, without duplicates
	foreach my $seengroup (sort (keys %seen)) {
		$failcount++ ;#####+=1;
# 		print "$RE has $failcount fails ($seengroup)\n<br>";	##### debugging
		
		my(@problems) = split (/\t/,$seen{$seengroup});
# 		print "$seengroup doesn't match with: @problems\n";	##### debugging
		foreach my $seenfail (@problems) { 
			my $AB = "$seengroup"."$seenfail";
			my $BA = "$seenfail"."$seengroup";
			next if ($duplicates{$AB});
			next if ($duplicates{$BA});
			
			if ($seengroup le $seenfail) {push @{$list{$RE}}, "$seengroup". "$seenfail";}
			else {push @{$list{$RE}}, "$seenfail". "$seengroup";}
						
# 			print "list(RE): @{$list{$RE}}\n";	##### debugging
			
			$duplicates{$AB} = 1;
			$duplicates{$BA} = 1;
		}
	}

		##### this checks for stringency, only enzymes that discriminate some fraction of groups can go on
		##### "list" has all enzymes, "partlist" has only passing enzymes
	if (defined $list{$RE}) {
# 		print "listRE; $list{$RE}\n<br>";
		@{$partlist{$RE}} = @{$list{$RE}} unless (scalar @{$list{$RE}} > $maxfail);
# 		print "partlist: @{$partlist{$RE}}\n" if (defined $partlist{$RE});	##### debugging
	}

	##### if there are any enzymes that singlehandedly distinguish all the groups, they are saved here
	if ($failcount == 0) {
		$success{$RE}=$site;
		@{$partlist{$RE}} = "NONE";	##### this will go to the yes_cuts files
		@{$list{$RE}} = "SUCCESS";	##### this will go ...?
	}
		
}

##### print to file the unlikely super-enzymes that singlehandedly distinguish all groups
if ($successfile) {
	open (SUCCESS, ">$successfile") or die "Couldn't open $successfile: $!";
	print SUCCESS "Enzyme\tcutsite\t (Success with fragments < $cutoff bp apart)";
	foreach my $successEnz (keys %success){
		print SUCCESS "\n$successEnz\t$success{$successEnz}";
	}
	close SUCCESS or die "Couldn't close $successfile: $!";
}
print ("<br>" . scalar keys %success);
print(" singly successful enzymes printed: <a href=\"./$timestring/success.txt\">success.txt</a>: enzymes that distinguish all the groups.\n<br>");


##### print to file the loser enzymes that didn't cut or didn't distinguish at all
if ($failfile) {
open (FAILURE, ">$failfile") or die "Couldn't open $failfile: $!";
print FAILURE "The following enzymes did not cut/distinguish in the region between $bpcutoff_low and $bpcutoff_high\n" . join "\n",(sort @nocutters);
close FAILURE or die "Couldn't open $failfile: $!";
}
print scalar @nocutters ." enzymes did not cut/distinguish in the region between $bpcutoff_low and $bpcutoff_high bp.  Print to <a href=\"./$timestring/no_cuts.txt\">no_cuts.txt</a>: enzymes that either do not cut or distinguish any groups\n<br>";


##### print to file the enzymes that cut and distinguished
if ($yescutsfile) {
	open (YESCUTS, ">$yescutsfile") or die "Couldn't open $yescutsfile: $!";
	print YESCUTS "#NDG\tEnzyme\tNon-discriminated groups";

	$listcount = keys %list;
	$partcount = keys %partlist;
	foreach my $okEnz (sort (keys %partlist)) {
		my $misscount = @{$partlist{$okEnz}};
		print YESCUTS "\n$misscount\t$okEnz\t$enzymes{$okEnz}\t". "@{$partlist{$okEnz}}";
	}
close YESCUTS or die "Couldn't close $yescutsfile: $!";
}
print "$partcount enzymes differentiate >" . (100*$stringency). "\% (".($numcombos-$maxfail) . ") group combinations, " . eval{$listcount - $partcount} . " enzymes do not.\n<br>";

##### print fragment lengths to file
if ($fragfile) {
	open (FRAGS, ">$fragfile") or die "Couldn't open $fragfile: $!";
	foreach my $taxaname (sort (keys %sequences)) {print FRAGS "\t$taxaname"}
	print FRAGS "\n";
	PRINTZYME: foreach my $printzyme (sort (keys %fragsave)) {
		
##		foreach my $allcutters (sort (keys %fragsave)) {	##### to print all fragments
##			next unless $allcutters =~ m/^$printzyme$/;
			print FRAGS "$printzyme\t" . (join "\t",@{$fragsave{$printzyme}}). "\n";
##			}
	}
close FRAGS or die "Couldn't close $fragfile: $!";
}
print "Fragments printed:<a href=\"./$timestring/fragfile.txt\">fragfile.txt</a>: fragment sizes for all enzymes\n<br>";











#######################
####### PART 2 ########  In which the qualities of each enzyme is compared with each other enzyme and
#######################  the best group of four enzymes comes out victorious.

my %zlist = matrixMaker (%partlist);
print "Enzyme match matrix printed:<a href=\"./$timestring/enzmatrix.txt\">enzmatrix.txt</a>: matrix of group combinations each passing enzyme cuts\n<br>";

# die "No enzyme groups will succeed because " . (join " and ",@zerogroup) . " are not discriminated by any enzymes\n" if (@zerogroup);
my %rev_groups = revMatrix (%zlist);
my $matchstring = 1 x $numcombos;
my $countprint = 1;
my $countprint2 = 1;
my $hashsize = keys(%rev_groups);
foreach my $dec_1 (sort {$a cmp $b} (keys %rev_groups)) {
		$hashsize--;
		if (($hashsize % 5) == 0) {print ("$hashsize")}
		else {print "."}
	foreach my $dec_2 (sort {$a cmp $b} (keys %rev_groups)) {
		next if $dec_2 le $dec_1;
		foreach my $dec_3 (sort {$a cmp $b} (keys %rev_groups)) {
			next if $dec_3 le $dec_2;
			foreach my $dec_4 (sort {$a cmp $b} (keys %rev_groups)) {
				next if $dec_4 le $dec_3;
				my $compare = ("$dec_1" | "$dec_2" | "$dec_3" | "$dec_4");
				my @goodsum = ($dec_1)=~ m/1/gi;    ### switch to tr// trick??
				push @goodsum,($dec_2)=~ m/1/gi;
				push @goodsum,($dec_3)=~ m/1/gi;
				push @goodsum,($dec_4) =~ m/1/gi;
				if ($compare eq $matchstring) {
					
					$combonumber++;
					push @{$goodcombo{$combonumber}},($rev_groups{$dec_1}[0],$rev_groups{$dec_2}[0],$rev_groups{$dec_3}[0],$rev_groups{$dec_4}[0]);
					$finalsave{$combonumber} = scalar @goodsum;
					
					
					
				}
			}
		}		
	}
}

if ($outfile) {
	open (OUT, ">$outfile") or die "Couldn't open $outfile: $!";
	print OUT "ENZYME PICKER KEY\nGroup Members\n\n";
	my %tetragroups;
	my $tetracount; 
	my $matchnumber;
	
	foreach my $printkey (keys %rev_groups) {
		$tetracount++;
		#####to print out isoschizimers
		foreach my $printkey3 (@{$rev_groups{$printkey}}) {
			push @{$tetragroups{$tetracount}},($printkey3, @{$isoschiz{$printkey3}});
		}
		
		print OUT "$tetracount\t@{$tetragroups{$tetracount}}\n";
		
	}

	##### a mess follows, but it works
	print OUT "\n\nENZYME SETS WITH FULL COVERAGE OF SEQUENCE GROUPS\nChoose from set above, below shows only first group member\nScore\tGroup A\tGroup B\tGroup C\tGroup D\n\n";
	
	##### cycle through all combos to order them by highest score
	foreach my $printkey2 (sort {$finalsave{$b} <=> $finalsave{$a}} (keys %finalsave)) {	##### keys are combonumbers, values are scalars
# 		print "printkey2: $printkey2\t$finalsave{$printkey2}\n";	##### debugging
		$matchnumber++;
		if ($matchnumber == $matchlimit) {last};
		printf OUT "\n%.2f", ($finalsave{$printkey2}/$numcombos);
		
		##### cycle through all four enzymes in goodcombo
		foreach my $printkey4 (sort {$a cmp $b} @{$goodcombo{$printkey2}}) {	##### cycling through 4 enzymes
			##### cycle through all the tetragroups in order to name them in finalout
			foreach my $printkey5 (keys %tetragroups) {		##### keys are tetracount numbers, values are enzymes+isoschizimers
				my @results = grep(/^$printkey4$/,@{$tetragroups{$printkey5}});
				if (@results) {print OUT "\t$printkey5 ($printkey4)\t";}
			}
			
		}
				
	}

	close OUT or die "Couldn't close $outfile: $!";
}

print "<br>$combonumber sets of enzymes were found that distinguish $numgroups groups.  The top $matchlimit were printed: <a href=\"./$timestring/finalout.txt\">finalout.txt</a>: list of 4-enzyme groups that will identify all groups\n<br>";

#####	___END MAIN___



# Print close of HTML file
print <<EOF ;
<h3>STOP</h3>
</body>
</html>
EOF

# <p>Your files are here:
# <br><a href="./$timestring/finalout.txt">finalout.txt</a>: list of 4-enzyme groups that will identify all groups
# <br><a href="./$timestring/success.txt">success.txt</a>: enzymes that distinguish all the groups
# <br><a href="./$timestring/no_cuts.txt">no_cuts.txt</a>: enzymes that do not cut/distinguish any groups
# <br><a href="./$timestring/enzmatrix.txt">enzmatrix.txt</a>: matrix of group combinations each passing enzyme cuts
# <br><a href="./$timestring/fragfile.txt">fragfile.txt</a>: fragment sizes for each passing enzyme
# <br><a href="./enzymes_comm_iso.txt">enzymes_comm_iso.txt</a>: the list of enzymes used for this analysis are all commercially available according to <a href="http://rebase.neb.com">REBASE</a>.


exit ;


















######################
#### SUBROUTINES #####
######################


# subroutine readEnzymeFile takes a tab delimited text file and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex

sub readEnzymeFile {
	my $inputFile = shift; 
	my %enzymes;
	open  (FH, $inputFile) or die "Couldn't open $inputFile: $!";
	while (my $line =<FH>) {
		$line =~ s/\r//;
		chomp $line;
		$line =~ m/ISO\t/;
		my ($RE, $site) = split (/\t/, $line);
		
		@{$isoschiz{$RE}} = split (/\t/,$');
		##### find out how many bases to add or subtract to get correct fragment lengths
		if ($site =~ s/\^//) {
			$cutlength{$RE} = length ($`);
			$enzymes{$RE} = uc ($site);
		}
		
		##### converts IUPAC code to regex
		$enzymes{$RE} =~ s/N/[ATCG]/g;
		$enzymes{$RE} =~ s/B/[TCG]/g;
		$enzymes{$RE} =~ s/D/[ATG]/g;
		$enzymes{$RE} =~ s/H/[ATC]/g;
		$enzymes{$RE} =~ s/V/[ACG]/g;
		$enzymes{$RE} =~ s/K/[GT]/g;
		$enzymes{$RE} =~ s/Y/[CT]/g;
		$enzymes{$RE} =~ s/S/[CG]/g;
		$enzymes{$RE} =~ s/W/[AT]/g;
		$enzymes{$RE} =~ s/M/[AC]/g;
		$enzymes{$RE} =~ s/R/[AG]/g;

	}
	
	close FH or die "Couldn't close $inputFile: $!";
	return %enzymes;
}

##### subroutine readFASTAfile takes a txt file in FASTA format and returns
##### a hash with each strain as a key and the uppercase sequence as a value

sub readFASTAfile {

#cgi 	my $inputFile = shift; 
	my $fastataxa;
	my %sequences;
	
	foreach (split("^", $cgivars{'fasta'})) {
		
		$line = sprintf("%s", $_);
# 		print "<br>line $line";
		$line =~ s/\r//;
		chomp $line;
		if ($line =~ m/^>/) {
		  	$fastataxa=$';
		}
		elsif ($line =~ m/^\s$/) {
			next;
		}
		else {
			$sequences{$fastataxa} .= uc ($line);
		}
		
	}

	return %sequences;
}

##### subroutine checkGroups [[checks for groups]]
#####

sub checkGroups {
# 	my $subset = ();
	foreach my $groupie (keys %sequences) {
		@taxasplit = split (/$splitter/,$groupie);
		my $grpid = @taxasplit[$splitnum];
# 		my $grpid = substr($groupie,$taxabid,$taxaeid);
# 		next if ($subset =~ m/($grpid)/);
		next if (grep(/^$grpid$/,@subs));
		push @subs, $grpid;
# 		$subset .= "$grpid\t";
	}
	
# 	@subs = split (/\t/, $subset);
	return @subs;
}

##### subroutine groupGroups makes array of all group combinations
#####
sub groupGroups {
	foreach my $sub1 (sort @subs) {
		foreach my $sub2 (sort @subs) {
			next if $sub2 le $sub1;
			push @grpgrp, $sub1.$sub2
		}
	}
# 	print "grpgrp: @grpgrp\n";	##### debugging
	return @grpgrp;
}

##### subroutine matrixMaker makes the 1/0 enzyme cut matrix
#####
sub matrixMaker {
open (MATRIX, ">$matrixfile") or die "Couldn't open $matrixfile: $!";
print MATRIX "\t" . (join "\t",@revgrpgrp) ."\n";
# print "grpgrp = @grpgrp\n<br>"; 
foreach my $enzgroup (sort (keys %partlist)) {
# 	print "MM: enzgroup: $enzgroup\n";	##### debugging
	push @{$zlist{$enzgroup}}, ("1") x $numcombos;
	my @zrow = @{$partlist{$enzgroup}};
# 	print "MM: zrow: @zrow\n";	##### debugging
		foreach my $zval (@zrow) {
# 			print "MM: zval: $zval\n<br>";	##### debugging
			next unless $zval;
			my $countup = 0;
			my $countdown = -1;
				foreach (@{$zlist{$enzgroup}}) {
					next unless defined;
					$zlist{$enzgroup}[$countdown] = "0" if $zval eq $grpgrp[$countup];
					$countup++;
					$countdown--;
				}
		}
	
	my(@matrixprint) = join "\t",@{$zlist{$enzgroup}};
	print MATRIX "$enzgroup\t@matrixprint\n";
	}
	
close MATRIX or die "Couldn't close $matrixfile: $!";
return %zlist;
}

##### subroutine revMatrix
##### has to do with increasing the efficiency by only running enzymes that have different cut patterns
#####
sub revMatrix {
	my %binary;
	foreach my $dbgroup (keys %zlist) {
# 		print "$dbgroup\n";	##### debugging
		$binary{$dbgroup} = join("",@{$zlist{$dbgroup}});
	}

		##### Reverse the hash to get just the unique values
		##### then scan through all the keys and add them to an array in the hash

	my %rev_dec = reverse %binary;
	foreach my $revkey (sort(keys %rev_dec)) {
		foreach my $forgroup (keys %binary) {
			if ($binary{$forgroup} eq $revkey) {
				push @{$rev_groups{$revkey}},$forgroup;
			}
		}
	}
	return %rev_groups;
}

#----------------- start of &getcgivars() module ----------------------
# borrowed from the internet
# Read all CGI vars into an associative array.
# If multiple input fields have the same name, they are concatenated into
#   one array element and delimited with the \0 character (which fails if
#   the input has any \0 characters, very unlikely but conceivably possible).
# Currently only supports Content-Type of application/x-www-form-urlencoded.
sub getcgivars {
    local($in, %in) ;
    local($name, $value) ;


    # First, read entire string of CGI vars into $in
    if ( ($ENV{'REQUEST_METHOD'} eq 'GET') ||
         ($ENV{'REQUEST_METHOD'} eq 'HEAD') ) {
        $in= $ENV{'QUERY_STRING'} ;

    } elsif ($ENV{'REQUEST_METHOD'} eq 'POST') {
        if ($ENV{'CONTENT_TYPE'}=~ m#^application/x-www-form-urlencoded$#i) {
            length($ENV{'CONTENT_LENGTH'})
                || &HTMLdie("No Content-Length sent with the POST request.") ;
            read(STDIN, $in, $ENV{'CONTENT_LENGTH'}) ;

        } else { 
            &HTMLdie("Unsupported Content-Type: $ENV{'CONTENT_TYPE'}") ;
        }

    } else {
        &HTMLdie("Script was called with unsupported REQUEST_METHOD.") ;
    }
    
    # Resolve and unencode name/value pairs into %in
    foreach (split(/[&;]/, $in)) {
#     	$OK_CHARS='-a-zA-Z0-9_.;|';	# A restrictive list, which
# 					# should be modified to match
# 					# an appropriate RFC, for example.
# 		s/[^$OK_CHARS]/_/go; ##insert for security

        s/\+/ /g ;
        ($name, $value)= split('=', $_, 2) ;
        $name=~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge ;
        $value=~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge ;
        $in{$name}.= "\0" if defined($in{$name}) ;  # concatenate multiple vars
        $in{$name}.= $value ;
    }

    return %in ;

}


# Die, outputting HTML error page
# If no $title, use a default title
sub HTMLdie {
    local($msg,$title)= @_ ;
    $title= "CGI Error" if $title eq '' ;
    print <<EOF ;
Content-type: text/html

<html>
<head>
<title>$title</title>
</head>
<body>
<h1>$title</h1>
<h3>$msg</h3>
</body>
</html>
EOF

    exit ;
}
