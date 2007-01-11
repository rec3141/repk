#!/usr/local/bin/perl -Tw
##### program to calculate length of terminal restriction fragments
##### for many enzymes and compare results among groups
##### input is a multiple sequence alignment in fasta format,
##### get rid of whitespace in sequences names and start sequences at beginning of primer

### (REQ FEATURE: split cutoff into INTERgroup and INTRAgroup)
### (REQ FEATURE: warn if any group is discriminated by no enzymes)
### (REQ FEATURE: print column of # groups each enzyme differentiates

$|++;    #forces a flush after every write or print on the currently selected output channel

use strict;
# use diagnostics;
use File::Copy;    #a module necessary for copying files

# use Data::Dumper;

# use CGI::Carp qw(fatalsToBrowser carpout);
# open( LOG, ">>./cgilogfile.repk.txt" ) or die("Unable to open log: $!\n");
# carpout(LOG);

######################################################
#                    PREAMBLE                        #
######################################################

##### First, print the CGI header
printCGI();
##### Then get the CGI variables into a list of strings
my %cgivars = &getcgivars;

##### Set up Results directory
my ( $sec, $min, $hour, $day, $mon, $year, undef, undef, undef ) = localtime time;
$year += 1900;
my $timestring = "RESULTS/REPK_" . "$year$mon$day-$hour$min-$sec";
mkdir( $timestring, 0755 );

##### variables set in web form
my $enzymeFile       = 'enzymes_type2.txt';             # default enzyme file
my $enzyme_list      = $cgivars{'enzyme_list'};         # selected enzymes
my $enzymeFileCustom = $cgivars{'enzymeFileCustom'};    # custom enzymes
my $fastafile        = $cgivars{'fasta'};               # FASTA formatted sequences
my $splitter         = $cgivars{'splitter'};            # character on which to split FASTA name
my $splitnum         = $cgivars{'splitnum'} - 1;        # the substr with the groups
my $cutoff           = $cgivars{'cutoff'};              # Cutoff length
my $bpcutoff_low     = $cgivars{'bpcutoff_low'};        # shortest acceptable fragment
my $bpcutoff_high    = $cgivars{'bpcutoff_high'};       # longest acceptable fragment
my $stringency       = $cgivars{'stringency'} / 100;    # Stringency filter
my $mismatches       = $cgivars{'mismatches'};          # max allowable group mismatches
my $matchlimit       = $cgivars{'matchlimit'};          # max matches to print, <1000
my $rdpinput         = $cgivars{'rdpinput'};            # output file from RDP classifier
( $matchlimit = 1000 ) if ( $matchlimit > 1000 );

##### warnings
print "<br><font color=red>WARNING</font>: <em>NO FASTA FILE DETECTED</em>" and exit unless $fastafile;
print "<br><font color=red>WARNING</font>: <em>CUTOFF</em> can not be greater than or equal to <em>BPCUTOFF_LOW</em>" and exit if ( $bpcutoff_low <= $cutoff );

##### Which enzymes were chosen
my $seenenz = 0;
my @enzymesList;
$seenenz++ if length($enzymeFileCustom) != 0;

if ( length($enzyme_list) != 0 ) {
    if ( $enzyme_list =~ m/allrebase/ ) { push( @enzymesList, 'allrebase' ); $seenenz++ }
    elsif ( $enzyme_list =~ m/norebase/ ) { push( @enzymesList, 'norebase' ) }
    else { @enzymesList = split( "\0", $enzyme_list ); $seenenz++ }
}

print "<br><font color=red>WARNING</font>: <em>NO ENZYMES CHOSEN</em>" and exit if ( $seenenz == 0 );

##### output filenames
my $renamedfile    = "./$timestring/renamed.fasta";     # renamed by RDP classifier
my $customfile     = "./$timestring/custom.txt";        # custom enzymes used or chosen
my $fragfile       = "./$timestring/fragfile.csv";      # fragment sizes for each passing enzyme
my $successfile    = "./$timestring/success.txt";       # enzymes that distinguish all the groups
my $yescutsfile    = "./$timestring/yes_cuts.txt";      # enzymes that cut and meet stringency
my $failfile       = "./$timestring/no_cuts.txt";       # enzymes that do not cut/distinguish any groups
my $matrixfile     = "./$timestring/enzmatrix.csv";     # matrix of group combinations each passing enzyme cuts
my $missingoutfile = "./$timestring/missingout.txt";    # partially successful combinations
my $outfile        = "./$timestring/finalout.txt";      # list of 4-enzyme groups that will identify all groups

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
my @justEnz;
my $combonumber     = 0;
my $usercombonumber = 0;
my $usermatchlimit  = $matchlimit;


##### Set up Sequences
my @renamed   = readRDPCfile( ( split( /\n/, $rdpinput ) ) );    # parse RDP-classifier output
my %sequences = readFASTAfile();                                 # parse FASTA sequences

my @subs = checkGroups(%sequences);                              # make the groups

my $numgroups = scalar(@subs);                                   # how many groups
my $numcombos = $numgroups * ( $numgroups - 1 ) / 2;             # how many combos of those groups
my $maxfail   = $numcombos - int( $numcombos * $stringency );    # how many to pass stringency

my @grpgrp    = groupGroups(@subs);                              # ???
my @revgrpgrp = reverse @grpgrp;                                 # ???

##### Print progress on Sequences
if ( $rdpinput =~ m/Classifier/g ) {
    open( RENAMED, ">$renamedfile" ) or die "Couldn't open $renamedfile: $!";
    foreach my $seqkey ( sort { $a cmp $b } keys %sequences ) {
        print RENAMED ">$seqkey\n$sequences{$seqkey}\n";
    }
    close RENAMED or die "Couldn't close $renamedfile: $!";
    print "PRINTED: Input sequences renamed by RDP-Classifier, <a href=\"$renamedfile\">renamed.fasta</a>\n<br>";
}

print "PROCESSING: There are $numgroups groups and $numcombos group combinations\n<br>";
print "PROCESSING: The minimum allowable group discriminations is $maxfail\n<br>";

##### Set up Enzymes
# If 'all enzymes' is chosen
if ( $enzymesList[0] eq 'allrebase' ) {
    %enzymes = readEnzfile($enzymeFile);
    print "PROCESSING: " . ( scalar keys %enzymes ) . " enzymes were read from <a href=\"$enzymeFile\">$enzymeFile</a>\n<br>";
}

# If 'all enzymes' is not chosen but others are
elsif ( $enzymesList[1] ) {
    %enzymes = readEnzchosen(@enzymesList);
    foreach (@enzymesList) { my ( $just, $junk ) = split( '\t', $_ ); push( @justEnz, $just ); }
    print "PROCESSING: " . scalar(@justEnz) . " enzymes were chosen from the list <a href=\"$customfile\">custom.txt</a>\n<br>";
}

if ($enzymeFileCustom) {
    %customenzymes = readEnzcustom($enzymeFileCustom);
    print "PROCESSING: " . ( scalar keys %customenzymes ) . " custom enzymes were used and added to <a href=\"$customfile\">custom.txt</a>\n<br>";
    @enzymes{ keys %customenzymes } = values %customenzymes;    # get

}

##### print custom enzymes
open( CUSTOM, ">$customfile" ) or die "Couldn't open $customfile: $!";
print CUSTOM join( "\n", ( @justEnz, ( keys %customenzymes ) ) ) . "\n";
close CUSTOM or die "Couldn't close $customfile: $!";

#######################
####### PART 1 ########  In which the fragments are cut with the enzymes, the terminal fragment lengths are saved,
#######################  the groups are differentiated, evil is vanquished, and peace is brought to the entire galaxy.

print "PROCESSING: Enzymes left to test ";
my $enzhashsize = keys(%enzymes);
ENZYME: foreach my $RE ( sort ( keys %enzymes ) ) {

    $enzhashsize--;
    if ( ( $enzhashsize % 10 ) == 0 ) { print("$enzhashsize") }
    else { print "." }

    ##### reverse complements each sequence to search for opposite strand cutsites
    $reverseEnzymes{$RE} = reverse( $enzymes{$RE} );
    $reverseEnzymes{$RE} =~ tr/ACGT[]/TGCA][/;

    my $site    = $enzymes{$RE};
    my $revsite = $reverseEnzymes{$RE};

    my %seen       = ();
    my %duplicates = ();
    my $failcount  = 0;

    ##### finds fragment lengths for all seqs, sets to negative value if not within bpcutoffs
    foreach my $taxa ( sort ( keys %sequences ) ) {
        $fragments{$taxa} = ();
        my $terminalFragment = 0;
        my $tempfrag         = 0;
        my $forwardFragment;
        my $reverseFragment;

        ##### checks the sequence for the forward recognition sequence, saves it to terminalFragment if passes stringency, tempfrag if not
        if ( $sequences{$taxa} =~ m/($site){1}/i ) {
            $forwardFragment = length($`) + $cutlength{$RE};
            if ( ( $forwardFragment < $bpcutoff_high ) and ( $forwardFragment > $bpcutoff_low ) ) {
                $terminalFragment = $forwardFragment;
            }
            else {
                $tempfrag = $forwardFragment;
            }

        }
        ##### checks the sequence for the reverse recognition sequence, saves it to terminalFragment if passes stringency
        ##### and is shorter than forward terminalFragment, tempfrag if it is shorter but doesn't pass stringency
        if ( $sequences{$taxa} =~ m/($revsite){1}/i ) {
            $reverseFragment = length($`) + $cutlength{$RE};
            if ( ( $reverseFragment < $bpcutoff_high ) and ( $reverseFragment > $bpcutoff_low ) and ( $reverseFragment < $forwardFragment ) ) {
                $terminalFragment = $reverseFragment;
            }
            else {
                ( $tempfrag = $reverseFragment ) if ( $reverseFragment < $tempfrag );
            }
        }

        ##### saves the fragment length that passes stringency for later printing and analysis
        if ($terminalFragment) {
            $fragments{$taxa} = "$terminalFragment";
            push @{ $fragsave{$RE} }, "$terminalFragment";
        }
        ##### saves the fragment length that does not pass stringency for printing later
        elsif ($tempfrag) {
            $fragments{$taxa} = 0;
            push @{ $fragsave{$RE} }, 1 * $tempfrag;
        }
        ##### saves the entire fragment length because there were no cuts
        else {
            $fragments{$taxa} = 0;
            push @{ $fragsave{$RE} }, 1 * length( $sequences{$taxa} );
        }

    }

    ##### saves any enzymes that cuts no sequences or cuts them all alike for later print out to failure file
    my %flip_frag = reverse %fragments;
    unless ( keys %flip_frag > 1 ) {
        push @nocutters, $RE;
        next ENZYME;
    }

    ##### compares groups, finds which groups have similar patterns
  CHECK: foreach my $checkkey ( sort ( keys %fragments ) ) {    ##### first keys are sequence names

        my @set_1array = split( /$splitter/, $checkkey );
        my $set_1      = $set_1array[$splitnum];

        foreach my $fragtaxa ( sort ( keys %fragments ) ) {     ##### second keys are sequence names
            next if ( $fragtaxa le $checkkey );                 ##### don't check the same taxa against itself or any less than it

            my @set_2array = split( /$splitter/, $fragtaxa );
            my $set_2 = $set_2array[$splitnum];
            next if ( $seen{$set_1} =~ m/($set_2)/ );           ##### don't check again if it's already been flagged
            my $lengthDifference = abs( $fragments{$checkkey} - $fragments{$fragtaxa} );    ##### find the difference in length between fragments

            next if ( $set_1 eq $set_2 );                                                   # we don't care as long as they're not different
            next if ( $lengthDifference > $cutoff );                                        # we DO want different groups to have different lengths

            $seen{$set_1} .= "$set_2\t";                                                    # flag the groups to know they won't match

        }
    }

    ##### saves which sets failed, without duplicates
    foreach my $seengroup ( sort ( keys %seen ) ) {
        $failcount++;

        my (@problems) = split( /\t/, $seen{$seengroup} );

        foreach my $seenfail (@problems) {
            my $AB = "$seengroup" . "$seenfail";
            my $BA = "$seenfail" . "$seengroup";
            next if ( $duplicates{$AB} );
            next if ( $duplicates{$BA} );

            if ( $seengroup le $seenfail ) { push @{ $list{$RE} }, "$seengroup" . "$seenfail"; }
            else { push @{ $list{$RE} }, "$seenfail" . "$seengroup"; }

            $duplicates{$AB} = 1;
            $duplicates{$BA} = 1;
        }
    }

    ##### this checks for stringency, only enzymes that discriminate some fraction of groups can go on
    ##### "list" has all enzymes, "partlist" has only passing enzymes
    if ( defined $list{$RE} ) {
        @{ $partlist{$RE} } = @{ $list{$RE} } unless ( scalar @{ $list{$RE} } > $maxfail );
    }

    ##### if there are any enzymes that singlehandedly distinguish all the groups, they are saved here
    if ( $failcount == 0 ) {
        $success{$RE} = $site;
        @{ $partlist{$RE} } = "NONE";       ##### this will go to the yes_cuts files
        @{ $list{$RE} }     = "SUCCESS";    ##### this will go ...?
    }

}

##############################
##### PRINT REULTS TO FILE   #
##############################

##### print to file the enzymes that cut and distinguished
if ($yescutsfile) {
    open( YESCUTS, ">$yescutsfile" ) or die "Couldn't open $yescutsfile: $!";
    print YESCUTS "#NDG\tEnzyme\tNon-discriminated groups";

    $listcount = keys %list;
    $partcount = keys %partlist;
    foreach my $okEnz ( sort ( keys %partlist ) ) {
        my $misscount = @{ $partlist{$okEnz} };
        print YESCUTS "\n$misscount\t$okEnz\t$enzymes{$okEnz}\t" . "@{$partlist{$okEnz}}";
    }
    close YESCUTS or die "Couldn't close $yescutsfile: $!";
}
print "<br>RESULTS: There were $partcount enzymes that differentiated >" . ( 100 * $stringency ) . "\% (" . ( $numcombos - $maxfail ) . ") of the group combinations, " . eval { $listcount - $partcount } . " enzymes did not.\n";

##### print fragment lengths to file
if ($fragfile) {
    open( FRAGS, ">$fragfile" ) or die "Couldn't open $fragfile: $!";
    foreach my $taxaname ( sort ( keys %sequences ) ) { print FRAGS "\t$taxaname" }
    print FRAGS "\n";
  PRINTZYME: foreach my $printzyme ( sort ( keys %fragsave ) ) {
        print FRAGS "$printzyme\t" . ( join "\t", @{ $fragsave{$printzyme} } ) . "\n";
    }
    close FRAGS or die "Couldn't close $fragfile: $!";
}
print "<br>PRINTED: Terminal restriction fragment lengths with all enzymes, <a href=\"./$timestring/fragfile.csv\">fragfile.csv</a>\n";

##### print to file the loser enzymes that didn't cut or didn't distinguish at all
if ($failfile) {
    open( FAILURE, ">$failfile" ) or die "Couldn't open $failfile: $!";
    print FAILURE "The following enzymes did not cut/distinguish in the region between $bpcutoff_low and $bpcutoff_high\n" . join "\n", ( sort @nocutters );
    close FAILURE or die "Couldn't open $failfile: $!";
}
print "<br>RESULTS: There were " . scalar @nocutters . " enzymes that did not cut or distinguish any groups in the region between $bpcutoff_low and $bpcutoff_high bp.\n";
print "<br>PRINTED: <a href=\"./$timestring/no_cuts.txt\">no_cuts.txt</a>\n" if (@nocutters);

##### print to file the unlikely super-enzymes that singlehandedly distinguish all groups
if ($successfile) {
    open( SUCCESS, ">$successfile" ) or die "Couldn't open $successfile: $!";
    print SUCCESS "Enzyme\tcutsite\t (Success with fragments < $cutoff bp apart)";
    foreach my $successEnz ( keys %success ) {
        print SUCCESS "\n$successEnz\t$success{$successEnz}";
    }
    close SUCCESS or die "Couldn't close $successfile: $!";
}
print "\n<br><em>RESULTS: There were " . ( scalar keys %success ) . " singly successful enzymes (enzymes that distinguish all the groups).\n</em><br>";

if ( scalar keys %success ) {
    print "PRINTED: <a href=\"./$timestring/success.txt\">success.txt</a>\n<br>";

    HTMLclose();
}

#######################
####### PART 2 ########  In which the qualities of each enzyme is compared with each other enzyme and
#######################  the best group of four enzymes comes out victorious.

my %zlist = matrixMaker(%partlist);

print "PRINTED: Enzyme match matrix, <a href=\"./$timestring/enzmatrix.csv\">enzmatrix.csv</a>\n<br>";
print "PROCESSING: Enzyme groups left to test ";

my %rev_groups = revMatrix(%zlist);

my $countprint  = 1;
my $countprint2 = 1;
my $hashsize    = keys(%rev_groups);

my @keylist;
for my $k1 ( sort { $a cmp $b } keys %rev_groups ) {
    push( @keylist, $k1 );
}

KEYLOOP: for my $i ( 0 .. $#keylist - 3 ) {
    $hashsize--;
    ( ( $hashsize % 5 ) == 0 ) ? print("$hashsize") : print(".");
    for my $j ( $i + 1 .. $#keylist - 2 ) {
        for my $k ( $j + 1 .. $#keylist - 1 ) {
            for my $l ( $k + 1 .. $#keylist ) {
                my $newcompare = ( "$keylist[$i]" | "$keylist[$j]" | "$keylist[$k]" | "$keylist[$l]" );
                my $newcount   = ( $newcompare =~ tr/1// );

                if ( $newcount == $numcombos ) {
                    $combonumber++;
                    if ( $combonumber > 9999 ) { print "<br><br><font color=red>WARNING</font>: <em><b>Whoa there!</b></em><br>Too many results, go back and try raising the enzyme stringency.<br><br>" and HTMLclose() }

                    push @{ $goodcombo{$combonumber} }, ( $rev_groups{ $keylist[$i] }[0], $rev_groups{ $keylist[$j] }[0], $rev_groups{ $keylist[$k] }[0], $rev_groups{ $keylist[$l] }[0] );
                    $finalsave{$combonumber} = ( ( $keylist[$i] ) =~ tr/1// ) + ( ( $keylist[$j] ) =~ tr/1// ) + ( ( $keylist[$k] ) =~ tr/1// ) + ( ( $keylist[$l] ) =~ tr/1// );
                }

                elsif ( $newcount >= ( $numcombos - $mismatches ) ) {
                    $usercombonumber++;
                    if ( $usercombonumber > 9999 ) { print "<br><br><font color=red>WARNING</font>: <em><b>Whoa there!</b></em><br>Too many results, go back and try lowering the number of allowable mismatches.<br><br>" and HTMLclose() }
                    push @{ $usercombo{$usercombonumber} }, ( $rev_groups{ $keylist[$i] }[0], $rev_groups{ $keylist[$j] }[0], $rev_groups{ $keylist[$k] }[0], $rev_groups{ $keylist[$l] }[0] );
                    $userfinalsave{$usercombonumber} = ( ( $keylist[$i] ) =~ tr/1// ) + ( ( $keylist[$j] ) =~ tr/1// ) + ( ( $keylist[$k] ) =~ tr/1// ) + ( ( $keylist[$l] ) =~ tr/1// );
                }
            }
        }
    }
}

if ($outfile) {
    open( OUT, ">$outfile" ) or die "Couldn't open $outfile: $!";
    print OUT "ENZYME PICKER KEY\nGroup Members\n\n";
    my %tetragroups;
    my $matchnumber = 0;
    my $tetracount;
    foreach my $printkey ( keys %rev_groups ) {
        $tetracount++;
        #####to print out isoschizimers
        foreach my $printkey3 ( @{ $rev_groups{$printkey} } ) {
            push @{ $tetragroups{$tetracount} }, ( $printkey3, @{ $isoschiz{$printkey3} } );
        }

        print OUT "$tetracount\t@{$tetragroups{$tetracount}}\n";

    }

    ##### a mess follows, but it works
    print OUT "\n\nENZYME SETS WITH FULL COVERAGE OF SEQUENCE GROUPS\nChoose from set above, below shows only first group member\nScore\tGroup A\tGroup B\tGroup C\tGroup D\n\n";

    ##### cycle through all combos to order them by highest score
    foreach my $printkey2 ( sort { $finalsave{$b} <=> $finalsave{$a} } ( keys %finalsave ) ) {    ##### keys are combonumbers, values are scalars

        if ( $matchnumber > $matchlimit ) { last }
        $matchnumber++;

        printf OUT "\n%.2f", ( $finalsave{$printkey2} / $numcombos );

        ##### cycle through all four enzymes in goodcombo
        foreach my $printkey4 ( sort { $a cmp $b } @{ $goodcombo{$printkey2} } ) {                ##### cycling through 4 enzymes
            ##### cycle through all the tetragroups in order to name them in finalout
            foreach my $printkey5 ( keys %tetragroups ) {                                         ##### keys are tetracount numbers, values are enzymes+isoschizimers
                my @results = grep( /^$printkey4$/, @{ $tetragroups{$printkey5} } );
                if (@results) { print OUT "\t$printkey5 ($printkey4)\t"; }
            }

        }

    }

    close OUT or die "Couldn't close $outfile: $!";
}

if ( $combonumber < 100 ) { $matchlimit = $combonumber }
if ( $combonumber > 0 ) {
    print "<br>RESULTS: There were $combonumber sets of enzymes found that distinguished $numgroups groups.\n";
    print "<br>PRINTED: The top $matchlimit combinations, <a href=\"./$timestring/finalout.txt\">finalout.txt</a>\n";
}
else { print "<br>RESULTS: <em>There were NO SUCCESSFUL ENZYME GROUPS, please try again with different parameters."; }



if ($missingoutfile) {
    open( USEROUT, ">$missingoutfile" ) or die "Couldn't open $missingoutfile: $!";
    print USEROUT "ENZYME PICKER KEY\nGroup Members\n\n";
    my %tetragroups;
    my $matchnumber = 0;
    my $tetracount;
    foreach my $printkey ( keys %rev_groups ) {
        $tetracount++;
        #####to print out isoschizimers
        foreach my $printkey3 ( @{ $rev_groups{$printkey} } ) {
            push @{ $tetragroups{$tetracount} }, ( $printkey3, @{ $isoschiz{$printkey3} } );
        }

        print USEROUT "$tetracount\t@{$tetragroups{$tetracount}}\n";

    }
    ##### a mess follows, but it works
    print USEROUT "\n\nENZYME SETS WITH FULL COVERAGE OF SEQUENCE GROUPS\nChoose from set above, below shows only first group member\nScore\tGroup A\tGroup B\tGroup C\tGroup D\n\n";

    ##### cycle through all combos to order them by highest score
    foreach my $printkey2 ( sort { $userfinalsave{$b} <=> $userfinalsave{$a} } ( keys %userfinalsave ) ) {    ##### keys are combonumbers, values are scalars

        $matchnumber++;
        if ( $matchnumber > $usermatchlimit ) { last }

        printf USEROUT "\n%.2f", ( $userfinalsave{$printkey2} / $numcombos );

        ##### cycle through all four enzymes in goodcombo
        foreach my $printkey4 ( sort { $a cmp $b } @{ $usercombo{$printkey2} } ) {                            ##### cycling through 4 enzymes
            ##### cycle through all the tetragroups in order to name them in finalout
            foreach my $printkey5 ( keys %tetragroups ) {                                                     ##### keys are tetracount numbers, values are enzymes+isoschizimers
                my @results = grep( /^$printkey4$/, @{ $tetragroups{$printkey5} } );
                if (@results) { print USEROUT "\t$printkey5 ($printkey4)\t"; }
            }

        }

    }

    close USEROUT or die "Couldn't close $outfile: $!";
}

if ( $usercombonumber < 100 ) { $usermatchlimit = $usercombonumber }
if ($mismatches) {
    if ( $usercombonumber > 0 ) {
        print "<br>RESULTS: There were $usercombonumber sets of enzymes found that distinguished at least " . ( $numgroups - $mismatches ) . " groups.\n";
        print "<br>PRINTED: The top $usermatchlimit partially successful combinations, <a href=\"./$timestring/missingout.txt\">missingout.txt</a>\n";
    }
    else { print "<br>RESULTS: <em>There were NO PARTIALLY SUCCESSFUL ENZYME GROUPS." }
}

#####	___END MAIN___

##### Print HTML close
HTMLclose();



######################
#### SUBROUTINES #####
######################

sub printCGI {

    # Print the CGI response header, required for all HTML output
    # Note the extra \n, to send the blank line
    print "Content-type: text/html\n\n";

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
}

# subroutine readEnzymeFile takes a tab delimited text file and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex
sub readEnzfile {
    my $inputFile = shift;

    open( FH, $inputFile ) or die "Couldn't open $inputFile: $!";
    my (@lines) = <FH>;    # read file into list
    close FH or die "Couldn't close $inputFile: $!";
    my %fileHash = readEnzymeText(@lines);
    return %fileHash;
}

sub readEnzcustom {
    my $line       = shift;
    my @lines      = split( /\n/, $line );
    my %customHash = readEnzymeText(@lines);
    return %customHash;
}

sub readEnzchosen {
    my %chosenAllHash = readEnzymeText(@_);
}

# subroutine readEnzymeCustom takes a CGI textarea input and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex

sub readEnzymeText {
    my %enzymeSub;
    foreach (@_) {

        s/\r//;
        chomp;
        my @REarray = split(/\s/);
        my $RE      = shift(@REarray);
        my $site    = shift(@REarray);
        shift(@REarray);    #kill the ISO
        @{ $isoschiz{$RE} } = @REarray;

        ##### finds out how many bases to add or subtract to get correct fragment lengths
        if ( $site =~ s/\^// ) {
            $cutlength{$RE} = length($`);
            $enzymeSub{$RE} = uc($site);
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
    foreach (@_) {
        if (m/Root;/) {
            ( my $sequence, my $results ) = split(/;Root;[0-9]{2,3};/);    # ';Root;100;'
            while ( $sequence =~ s/\W;[;]$// ) { chop($sequence) }
            chop $sequence;
            my @taxa = split( /;/, $results );
            grep( s/\s*// =~ $_, @taxa );
            push( @renarray, join( '_', ( @taxa[ 0, 2, 4, 6, 8, 10 ], $sequence ) ) );
        }
    }
    return @renarray;
}

##### subroutine readFASTAfile takes a txt file in FASTA format and returns
##### a hash with each strain as a key and the uppercase sequence as a value

sub readFASTAfile {

    my $fastataxa;
    my %sequences;

    foreach ( split( "^", $fastafile ) ) {
        my $line = sprintf( "%s", $_ );
        $line =~ s/\r//;    # should be taken care of by CGI security
        chomp $line;

        if ( $line =~ m/^\s$/ ) { next }
        elsif ( $line =~ m/FASTASTART/ ) { (@renamed) ? ( $fastataxa = shift(@renamed) ) : ( $fastataxa = $' ) }
        else {
            $line =~ s/[-.]//g;    # remove gaps from sequences
            $sequences{$fastataxa} .= uc($line);    # save
        }
    }
    return %sequences;
}

##### subroutine checkGroups [[checks for groups]]
#####

sub checkGroups {
    foreach my $groupie ( keys %sequences ) {
        my @taxasplit = split( /$splitter/, $groupie );
        print( "<br><font color=red>WARNING</font>: Not enough groups -- lower your group subset to " . scalar(@taxasplit) . " or less" ) and exit if ( scalar(@taxasplit) <= $splitnum );
        my $grpid = $taxasplit[$splitnum];

        next if ( grep( /^$grpid$/, @subs ) );
        push @subs, $grpid;
    }
    return @subs;
}

##### subroutine groupGroups makes array of all group combinations
#####
sub groupGroups {
    foreach my $sub1 ( sort @subs ) {
        foreach my $sub2 ( sort @subs ) {
            next if $sub2 le $sub1;
            push @grpgrp, $sub1 . $sub2;
        }
    }

    return @grpgrp;
}

##### subroutine matrixMaker makes the 1/0 enzyme cut matrix
#####

sub matrixMaker {
    open( MATRIX, ">$matrixfile" ) or die "Couldn't open $matrixfile: $!";
    print MATRIX "\t" . ( join "\t", @revgrpgrp ) . "\n";

    foreach my $enzgroup ( sort ( keys %partlist ) ) {

        push @{ $zlist{$enzgroup} }, ("1") x $numcombos;
        my @zrow = @{ $partlist{$enzgroup} };

        foreach my $zval (@zrow) {
            next unless $zval;
            my $countup   = 0;
            my $countdown = -1;
            foreach ( @{ $zlist{$enzgroup} } ) {
                next unless defined;
                $zlist{$enzgroup}[$countdown] = "0" if $zval eq $grpgrp[$countup];
                $countup++;
                $countdown--;
            }
        }

        my (@matrixprint) = join "\t", @{ $zlist{$enzgroup} };
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
    foreach my $dbgroup ( keys %zlist ) {

        $binary{$dbgroup} = join( "", @{ $zlist{$dbgroup} } );
    }

    ##### Reverse the hash to get just the unique values
    ##### then scan through all the keys and add them to an array in the hash

    my %rev_dec = reverse %binary;
    foreach my $revkey ( sort( keys %rev_dec ) ) {
        foreach my $forgroup ( keys %binary ) {
            if ( $binary{$forgroup} eq $revkey ) {
                push @{ $rev_groups{$revkey} }, $forgroup;
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
# Print the CGI variables sent by the user.
# Note that the order of variables is unpredictable.
# Also note this simple example assumes all input fields had unique names,
#   though the &getcgivars() routine correctly handles similarly named
#   fields-- it delimits the multiple values with the \0 character, within
#   $cgivars{$_}.
# foreach (keys %cgivars) {
#    print "<li>[$_] = [$cgivars{$_}]\n" ;
# }

sub getcgivars {
    my $in,   my %in;
    my $name, my $value;

    # First, read entire string of CGI vars into $in
    if ( ( $ENV{'REQUEST_METHOD'} eq 'GET' ) || ( $ENV{'REQUEST_METHOD'} eq 'HEAD' ) ) { $in = $ENV{'QUERY_STRING'} }
    elsif ( $ENV{'REQUEST_METHOD'} eq 'POST' ) {
        if ( $ENV{'CONTENT_TYPE'} =~ m#^application/x-www-form-urlencoded$#i ) {
            length( $ENV{'CONTENT_LENGTH'} ) || &HTMLdie("No Content-Length sent with the POST request.");
            read( STDIN, $in, $ENV{'CONTENT_LENGTH'} );
        }
        else { &HTMLdie("Unsupported Content-Type: $ENV{'CONTENT_TYPE'}") }
    }
    else { &HTMLdie("Script was called with unsupported REQUEST_METHOD.") }

    # Resolve and unencode name/value pairs into %in
    foreach ( split( /[&;]/, $in ) ) {

        s/\+/ /g;    # ???
        ( $name, $value ) = split( '=', $_, 2 );
        $name  =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge;    # convert names to proper characters
        $value =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge;    # convert values to proper characters


        # CGI security, removes dangerous characters
        if ($name !~ m/fasta/) {$value =~ s/[^-a-zA-Z0-9_.;|\n]//g}
        elsif ($name =~ m/fasta/) {$value =~ s/^\>|\n\>/^FASTASTART/g;$value =~ s/[^-a-zA-Z0-9_.;|\n]//g}

        $in{$name} .= "\0" if defined( $in{$name} );      # concatenate multiple vars
        $in{$name} .= $value;
    }

    return %in;

}

# Die, outputting HTML error page
# If no $title, use a default title
sub HTMLdie {
    my $msg, my $title = @_;
    $title = "CGI Error" if $title eq '';
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

    exit;
}

sub HTMLclose {

    # Print close of HTML file
    print <<EOF ;
<h3>DONE</h3>
</body>
</html>
EOF

    exit;
}
