#!/usr/local/bin/perl -Tw
# use strict;
use diagnostics;
$|++; #forces a flush after every write or print on the currently selected output channel

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

use CGI::Carp qw(fatalsToBrowser carpout);
 open(LOG, ">>./cgilogfile.repk.txt") or
 die("Unable to open log: $!\n");
 carpout(LOG);

use File::Copy;	#a module necessary for copying files
use Data::Dumper;


######################
######################
#####################

my($sec,$min,$hour,$day,$mon,$year,undef,undef,undef) = localtime time;$year += 1900;
my $timestring = "RESULTS/REPK_"."$year$mon$day-$hour$min-$sec"; 
mkdir($timestring, 0755);


# First, get the CGI variables into a list of strings
my %cgivars= &getcgivars ;

# Print the CGI response header, required for all HTML output
# Note the extra \n, to send the blank line
print "Content-type: text/html\n\n" ;

# Finally, print out the complete HTML response page
print <<EOF ;
<html>
<head><title>Repicker</title></head>
<body>
<h2>Results Page</h2>
<br>Links can be opened into new tabs or saved to disk by right-clicking,
<br>but <em><b>do not click directly on the files</b></em> until the run is over.
<br><br>You can access these files for up to 48 hours by bookmarking <a href="http://staff.washington.edu/rec3141/repk/$timestring/">this URL</a><br><br>
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


##### variables that can be set
my $rdpinput = $cgivars{'rdpinput'}; #from RDP classifier
my $cutoff = $cgivars{'cutoff'}; #5;	##### furthest apart two fragments can be in length and still be considered the same fragments
	##### could split up the cutoff into 2: one for INTERgroup, another for INTRAgroup
my $stringency = $cgivars{'stringency'}/100; #0.0; 	# an enzyme must distinguish MORE than this percentage of groups to be acceptable
#my $matrixwarning = $cgivars{'matrixwarning'}; #0.1;# doesn't work, but should: display warning about groups that are cut by fewer than this fraction of enzymes
my $bpcutoff_low = $cgivars{'bpcutoff_low'}; #75;	# shortest acceptable fragment
my $bpcutoff_high = $cgivars{'bpcutoff_high'}; #900;# longest acceptable fragment
my $matchlimit = $cgivars{'matchlimit'}; #100;	# the most number of matches we want printed
($matchlimit = 1000) if ($matchlimit > 1000);
my $usermatchlimit = $matchlimit;

my $splitter = $cgivars{'splitter'}; #"_";
my $splitnum = $cgivars{'splitnum'} -1;
my $enzymeFileCustom = $cgivars{'enzymeFileCustom'} if (defined $cgivars{'enzymeFileCustom'}); #custom enzymes
if ((defined $cgivars{'enzyme_list'}) && ($cgivars{'enzyme_list'} ne 'norebase')) {

 my $enzymeChosen = $cgivars{'enzyme_list'} ; #selected enzymes
 @enzymesList = split("\0",$enzymeChosen);
grep{s/norebase//} @enzymesList;
# print join('-',@enzymesList);
}
else {print "<br>PROCESSING: <em>NO ENZYMES CHOSEN" and exit unless $enzymeFileCustom}


my $mismatches = $cgivars{'mismatches'};

##### input filenames
my $enzymeFile = 'enzymes_type2.txt';			#"EnzymeInputFile.txt";

##### this program makes copies of the enzyme file and sequence alignment file
#copy("$enzymeFile","./"."$timestring"."/");
#copy("$alignmentFile","./"."$timestring"."/");

# ##### output filenames
my $yescutsfile = "./$timestring/yes_cuts.txt";		# enzymes that cut and meet stringency
my $successfile = "./$timestring/success.txt";	# enzymes that distinguish all the groups
my $failfile = "./$timestring/no_cuts.txt";		# enzymes that do not cut/distinguish any groups
my $matrixfile = "./$timestring/enzmatrix.csv";	# matrix of group combinations each passing enzyme cuts
my $outfile = "./$timestring/finalout.txt";		# list of 4-enzyme groups that will identify all groups
my $fragfile = "./$timestring/fragfile.csv";		# fragment sizes for each passing enzyme
my $customfile = "./$timestring/custom.txt";	#custom enzymes used or chosen
my $renamedfile = "./$timestring/renamed.fasta"; #renamed by RDP classifier
my $missingoutfile = "./$timestring/missingout.txt"; #partially successful combinations

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
my %usercombo;
my %enzymes;
my %customenzymes;
my %userfinalsave;
my %finalsave;

##### variables that are printed out of scope
my $listcount;
my $partcount;
my @nocutters;
my $combonumber=0;
my $usercombonumber=0;
my @subs;
my @justEnz;
# my @zerogroup;

#######################
####### PART 1 ########  In which the fragments are cut with the enzymes, the terminal fragment lengths are saved,
#######################  the groups are differentiated, evil is vanquished, and peace is brought to the entire galaxy.

my @renamed = readRDPCfile((split(/\n/, $rdpinput)));
my %sequences = readFASTAfile(%cgivars); #cgi $alignmentFile);

if ($rdpinput =~ m/Classifier/g) {
 open (RENAMED, ">$renamedfile") or die "Couldn't open $renamedfile: $!";
 foreach my $seqkey (sort {$a cmp $b} keys %sequences) {
  print RENAMED ">$seqkey\n$sequences{$seqkey}\n";
 }
 close RENAMED or die "Couldn't close $renamedfile: $!";
 print "PRINTED: Input sequences renamed by RDP-Classifier, <a href=\"./$timestring/renamed.fasta\">renamed.fasta</a>\n<br>";
}

@subs = checkGroups (%sequences);

my $numgroups = @subs;
my $numcombos = $numgroups*($numgroups-1)/2;
print "PROCESSING: There are $numgroups groups and $numcombos group combinations\n<br>";
my $maxfail = $numcombos - int($numcombos*$stringency);
print "PROCESSING: The minimum allowable group discriminations is $maxfail\n<br>";

my @grpgrp = groupGroups (@subs);
my @revgrpgrp = reverse @grpgrp;




if (grep($_ =~ /allrebase/, @enzymesList)) {
%enzymes = itsafile($enzymeFile);
print "PROCESSING: ". (scalar keys %enzymes) . " enzymes were read from <a href=\"$enzymeFile\">$enzymeFile</a>\n<br>";
}

elsif (($_ !~ /allrebase/) && (defined @enzymesList))  {

%enzymes = itsachosen(@enzymesList);
# print @enzymesList;
foreach (@enzymesList) {my($just,$junk) = split('\t',$_); push(@justEnz,$just)}
print "PROCESSING: ". scalar(@enzymesList) . " enzymes were chosen from the list <a href=\"$customfile\">custom.txt</a>\n<br>";
}

if ($enzymeFileCustom) {
%customenzymes = itsacustom($enzymeFileCustom);
print "PROCESSING: ". (scalar keys %customenzymes) . " custom enzymes were used <a href=\"$customfile\">custom.txt</a>\n<br>";
@enzymes{keys %customenzymes} = values %customenzymes;
}

open (CUSTOM, ">$customfile") or die "Couldn't open $customfile: $!";
print CUSTOM join("\n",(@justEnz ,(keys %customenzymes))) ."\n";
close CUSTOM or die "Couldn't close $customfile: $!";


#   Groups distinguished by fewer than " . (scalar keys %enzymes) * $matrixwarning . " enzymes will be printed to screen.\n<br>");

print "PROCESSING: Enzymes left to test ";
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
		my $set_1 = $set_1array[$splitnum];
		
		foreach my $fragtaxa (sort (keys %fragments)) {	##### second keys are sequence names
			next if ($fragtaxa le $checkkey);			##### don't check the same taxa against itself or any less than it
# 			my $set_2 = substr ($fragtaxa, $taxabid, $taxaeid);	##### pulls group name out of second sequence name
			my @set_2array = split (/$splitter/, $fragtaxa);
			my $set_2 = $set_2array[$splitnum];
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
print "<br>RESULTS: There were $partcount enzymes that differentiated >" . (100*$stringency). "\% (".($numcombos-$maxfail) . ") of the group combinations, " . eval{$listcount - $partcount} . " enzymes did not.\n";

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
print "<br>PRINTED: Terminal restriction fragment lengths with all enzymes, <a href=\"./$timestring/fragfile.csv\">fragfile.csv</a>\n";


##### print to file the loser enzymes that didn't cut or didn't distinguish at all
if ($failfile) {
open (FAILURE, ">$failfile") or die "Couldn't open $failfile: $!";
print FAILURE "The following enzymes did not cut/distinguish in the region between $bpcutoff_low and $bpcutoff_high\n" . join "\n",(sort @nocutters);
close FAILURE or die "Couldn't open $failfile: $!";
}
print "<br>RESULTS: There were " . scalar @nocutters ." enzymes that did not cut or distinguish any groups in the region between $bpcutoff_low and $bpcutoff_high bp.\n";
print "<br>PRINTED: <a href=\"./$timestring/no_cuts.txt\">no_cuts.txt</a>\n" if (@nocutters);


##### print to file the unlikely super-enzymes that singlehandedly distinguish all groups
if ($successfile) {
	open (SUCCESS, ">$successfile") or die "Couldn't open $successfile: $!";
	print SUCCESS "Enzyme\tcutsite\t (Success with fragments < $cutoff bp apart)";
	foreach my $successEnz (keys %success){
		print SUCCESS "\n$successEnz\t$success{$successEnz}";
	}
	close SUCCESS or die "Couldn't close $successfile: $!";
}
print "\n<br><em>RESULTS: There were " . (scalar keys %success) . " singly successful enzymes (enzymes that distinguish all the groups).\n</em><br>";

if (scalar keys %success) {
 print "PRINTED: <a href=\"./$timestring/success.txt\">success.txt</a>\n<br>";

# Print close of HTML file
print <<EOF ;
<h3>DONE</h3>
</body>
</html>
EOF

exit ;
}




#######################
####### PART 2 ########  In which the qualities of each enzyme is compared with each other enzyme and
#######################  the best group of four enzymes comes out victorious.

my %zlist = matrixMaker (%partlist);

print "PRINTED: Enzyme match matrix, <a href=\"./$timestring/enzmatrix.csv\">enzmatrix.csv</a>\n<br>";
print "PROCESSING: Enzymes left to test ";
# die "No enzyme groups will succeed because " . (join " and ",@zerogroup) . " are not discriminated by any enzymes\n" if (@zerogroup);
my %rev_groups = revMatrix (%zlist);
# print Dumper %rev_groups;

my $countprint = 1;
my $countprint2 = 1;
my $hashsize = keys(%rev_groups);

my @keylist;
for my $k1 (sort {$a cmp $b} keys %rev_groups) {
push(@keylist,$k1);
# print "keys: " . join("\t",@keylist) . "\n";
}

KEYLOOP: for my $i (0..$#keylist-3) { 
 $hashsize--;
 (($hashsize % 5) == 0) ? print("$hashsize") : print(".");
 for my $j ($i+1..$#keylist-2) {
  for my $k ($j+1..$#keylist-1) {
   for my $l ($k+1..$#keylist) {
	my $newcompare = ("$keylist[$i]" | "$keylist[$j]" | "$keylist[$k]" | "$keylist[$l]");
	my $newcount = ($newcompare =~ tr/1//);



	if ($newcount == $numcombos) {
	 $combonumber++;
		if ($combonumber > 9999) {print "<br><br><em><b>Whoa there!</b></em><br>Too many results, go back and raise the enzyme stringency.  <br>Stopping soon...<br><br>" and last KEYLOOP}

	 push @{$goodcombo{$combonumber}},($rev_groups{$keylist[$i]}[0],$rev_groups{$keylist[$j]}[0],$rev_groups{$keylist[$k]}[0],$rev_groups{$keylist[$l]}[0]);
	 $finalsave{$combonumber} = (($keylist[$i])=~ tr/1//) + (($keylist[$j])=~ tr/1//) + (($keylist[$k])=~ tr/1//) + (($keylist[$l])=~ tr/1//);
	}

	elsif ($newcount >= ($numcombos - $mismatches)) {
	 $usercombonumber++;
		if ($usercombonumber > 9999) {print "<br><br><em><b>Whoa there!</b></em><br>Too many results, go back and lower the number of allowable mismatches.  <br>Stopping soon...<br><br>" and last KEYLOOP}
	 push @{$usercombo{$usercombonumber}},($rev_groups{$keylist[$i]}[0],$rev_groups{$keylist[$j]}[0],$rev_groups{$keylist[$k]}[0],$rev_groups{$keylist[$l]}[0]);
	 $userfinalsave{$usercombonumber} = (($keylist[$i])=~ tr/1//) + (($keylist[$j])=~ tr/1//) + (($keylist[$k])=~ tr/1//) + (($keylist[$l])=~ tr/1//);
	}
   }
  }
 }
}



#SHOULD PRINT A COLUMN OF HOW MANY GROUPS EACH ENZYME DIFFERENTIATED
if ($outfile) {
	open (OUT, ">$outfile") or die "Couldn't open $outfile: $!";
	print OUT "ENZYME PICKER KEY\nGroup Members\n\n";
	my %tetragroups;
	my $matchnumber=0;
	my $tetracount; 
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
		if ($matchnumber > $matchlimit) {last};
		$matchnumber++;



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

if ($combonumber < 100) {$matchlimit = $combonumber};
if ($combonumber > 0) {
	print "<br>RESULTS: There were $combonumber sets of enzymes found that distinguished $numgroups groups.\n";
	print "<br>PRINTED: The top $matchlimit combinations, <a href=\"./$timestring/finalout.txt\">finalout.txt</a>\n"; 
	}
else {print "<br>RESULTS: <em>There were NO SUCCESSFUL ENZYME GROUPS, please try again with different parameters.";}

# print "<br>RESULTS: There were $combonumber sets of enzymes found that distinguished $numgroups groups.\n<br>";
# print "PRINTED: The top $matchlimit combinations, <a href=\"./$timestring/finalout.txt\">finalout.txt</a>\n<br>" if ($combonumber > 0);


if ($missingoutfile) {
	open (USEROUT, ">$missingoutfile") or die "Couldn't open $missingoutfile: $!";
	print USEROUT "ENZYME PICKER KEY\nGroup Members\n\n";
	my %tetragroups;
	my $matchnumber=0;
	my $tetracount;
	foreach my $printkey (keys %rev_groups) {
		$tetracount++;
		#####to print out isoschizimers
		foreach my $printkey3 (@{$rev_groups{$printkey}}) {
			push @{$tetragroups{$tetracount}},($printkey3, @{$isoschiz{$printkey3}});
		}
		
		print USEROUT "$tetracount\t@{$tetragroups{$tetracount}}\n";
		
	}
	##### a mess follows, but it works
	print USEROUT "\n\nENZYME SETS WITH FULL COVERAGE OF SEQUENCE GROUPS\nChoose from set above, below shows only first group member\nScore\tGroup A\tGroup B\tGroup C\tGroup D\n\n";

	##### cycle through all combos to order them by highest score
	foreach my $printkey2 (sort {$userfinalsave{$b} <=> $userfinalsave{$a}} (keys %userfinalsave)) {	##### keys are combonumbers, values are scalars
# 		print "printkey2: $printkey2\t$finalsave{$printkey2}\n";	##### debugging
		$matchnumber++;
		if ($matchnumber > $usermatchlimit) {last};

		printf USEROUT "\n%.2f", ($userfinalsave{$printkey2}/$numcombos);
		
		##### cycle through all four enzymes in goodcombo
		foreach my $printkey4 (sort {$a cmp $b} @{$usercombo{$printkey2}}) {	##### cycling through 4 enzymes
			##### cycle through all the tetragroups in order to name them in finalout
			foreach my $printkey5 (keys %tetragroups) {		##### keys are tetracount numbers, values are enzymes+isoschizimers
				my @results = grep(/^$printkey4$/,@{$tetragroups{$printkey5}});
				if (@results) {print USEROUT "\t$printkey5 ($printkey4)\t";}
			}

		}
				
	}

	close USEROUT or die "Couldn't close $outfile: $!";
}

if ($usercombonumber < 100) {$usermatchlimit = $usercombonumber};
if ($mismatches) {
 if ($usercombonumber > 0) {
	print "<br>RESULTS: There were $usercombonumber sets of enzymes found that distinguished at least " . ($numgroups-$mismatches) . " groups.\n";
	print "<br>PRINTED: The top $usermatchlimit partially successful combinations, <a href=\"./$timestring/missingout.txt\">missingout.txt</a>\n";
	}
	else {print "<br>RESULTS: <em>There were NO PARTIALLY SUCCESSFUL ENZYME GROUPS."}
}





#####	___END MAIN___

# Print close of HTML file
print <<EOF ;
<h3>DONE</h3>
</body>
</html>
EOF

exit ;


















######################
#### SUBROUTINES #####
######################


# subroutine readEnzymeFile takes a tab delimited text file and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex
sub itsafile {
my $inputFile = shift;

open  (FH, $inputFile) or die "Couldn't open $inputFile: $!";
my(@lines) = <FH>; # read file into list
close FH or die "Couldn't close $inputFile: $!";
my %fileHash = readEnzymeText(@lines);
return %fileHash;
}

sub itsacustom {
my $line = shift;
my @lines = split(/\n/,$line);
my %customHash = readEnzymeText(@lines);
return %customHash
}

sub itsachosen {
my %chosenAllHash = readEnzymeText(@_);
# print Dumper %chosenAllHash;
}



# subroutine readEnzymeCustom takes a CGI textarea input and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex

sub readEnzymeText {
	my %enzymeSub;
	foreach (@_) {
# 	print "<br>loopy: $_";
		s/\r//;
		chomp;
		my ($RE, $site) = split(/\s/);   #changed from \t to \s
		
		@{$isoschiz{$RE}} = split (/\s/,$'); #changed from \t to \s
		##### finds out how many bases to add or subtract to get correct fragment lengths
		if ($site =~ s/\^//) {
			$cutlength{$RE} = length ($`);
			$enzymeSub{$RE} = uc ($site);
		}
		
		##### converts IUPAC code to regex
		$enzymeSub{$RE} =~ s/N/[ATCG]/g;
		$enzymeSub{$RE} =~ s/B/[TCG]/g;
		$enzymeSub{$RE} =~ s/D/[ATG]/g;
		$enzymeSub{$RE} =~ s/H/[ATC]/g;
		$enzymeSub{$RE} =~ s/V/[ACG]/g;
		$enzymeSub{$RE} =~ s/K/[GT]/g;
		$enzymeSub{$RE} =~ s/Y/[CT]/g;
		$enzymeSub{$RE} =~ s/S/[CG]/g;
		$enzymeSub{$RE} =~ s/W/[AT]/g;
		$enzymeSub{$RE} =~ s/M/[AC]/g;
		$enzymeSub{$RE} =~ s/R/[AG]/g;

}

	return %enzymeSub;
}

sub readRDPCfile {
my @renarray;
foreach(@_) {
  if (m/Root;/) {
   (my $sequence,my $results) = split(/Root;\s+[0-9]{2,3}[\%];/); # 'Root; 100%;'
   while ($sequence =~ s/\W;[;\s]$//) {chop ($sequence)}
   chop $sequence;
   my @taxa = split(/;/, $results);
   grep(s/\s*// =~ $_, @taxa);
   push(@renarray, join('_',(@taxa[0,2,4,6,8,10], $sequence)));
  }
 }
# print Dumper @renarray;
return @renarray;
}


##### subroutine readFASTAfile takes a txt file in FASTA format and returns
##### a hash with each strain as a key and the uppercase sequence as a value

sub readFASTAfile {

#cgi 	my $inputFile = shift; 
	my $fastataxa;
	my %sequences;
	
	foreach (split("^", $cgivars{'fasta'})) {
		
		my $line = sprintf("%s", $_);
# 		print "<br>line $line";
		$line =~ s/\r//;
		chomp $line;
		if ($line =~ m/^>/) { (@renamed) ? ($fastataxa=shift(@renamed)) : ($fastataxa=$');
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
		my @taxasplit = split (/$splitter/,$groupie);
		my $grpid = $taxasplit[$splitnum];
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
# print Dumper %zlist;
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
# by James Marshall (http://www.jmarshall.com/easy/cgi/hello.pl.txt)
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
