#!/opt/local/bin/perl
# -w
#-T
##### program to calculate length of terminal restriction fragments
##### for many enzymes and compare results among groups
##### input is a multiple sequence alignment in fasta format,
##### get rid of whitespace in sequences names and start sequences at beginning of primer

### (REQ FEATURE: split cutoff into INTERgroup and INTRAgroup)

$|++;    #forces a flush after every write or print on the currently selected output channel

use strict;
#use diagnostics;
use File::Copy;    #a module necessary for copying files


# use Data::Dumper;

use CGI::Carp qw(fatalsToBrowser carpout);
#open(LOG, ">>/Library/Webserver/Documents/cgi/repk/RESULTS/cgilogfile.repk.txt" ) or die("Unable to open log: $!\n");
#carpout(LOG);

######################################################
#                    PREAMBLE                        #
######################################################

##### Set up Results directory
my ( $sec, $min, $hour, $day, $mon, $year, undef, undef, undef ) = localtime time;
$year += 1900;
my $rand_hex = join "", map { unpack "H*", chr(rand(128)) } 1..8;
my $tempstring = "REPK_" . "$year$mon$day-$hour$min-$rand_hex";  #;" . int(rand(1000));
my $timestring = "RESULTS/" . $tempstring;
mkdir( $timestring, 0755 );

##### Set up file locations
my $perl_location = "/opt/local/bin/perl";
my $prettymatrix_location = "/usr/local/bin/prettymatrix3.pl";
my $repk_location = "/Library/Webserver/Documents/cgi/repk/";
my $cgi_location="/cgi/repk/";

#### write .htaccess to allow indexing ####
open(HTACCESS, ">>./$timestring/.htaccess") or die "couldn't open .htaccess: $!";
print HTACCESS "Options +Indexes";
close(HTACCESS);

##### Then get the CGI variables into a list of strings
my %cgivars = ();
%cgivars = &getcgivars;

##### First, print the CGI header
printCGI();

##### for repeating if too low stringency 
my $repeats = 0;

BEGINNING:
if ($repeats == 0) {open( INDEX, ">./$timestring/results.html" ) or die "Couldn't open results.html: $!";}
else {close INDEX; open( INDEX, ">./$timestring/results.html" ) or die "Couldn't open results.html: $!";}

printResults() unless ($cgivars{'taxacheck'});

# use GD::Simple;
# print `perl -v`;





##### variables set in web form
my $enzymeFile       = 'enzymes_type2p.txt';             # default enzyme file
my $enzyme_list      = $cgivars{'enzyme_list'};         # selected enzymes
my $enzymeFileCustom = $cgivars{'enzymeFileCustom'};    # custom enzymes
my $fastafile        = $cgivars{'fasta'};               # FASTA formatted sequences
# my $splitter         = $cgivars{'splitter'};            # character on which to split FASTA name
my $splitter         = "_";                             # character on which to split FASTA name
my $splitnum         = $cgivars{'splitnum'} - 1;        # the substr with the groups
my $taxacheck        = $cgivars{'taxacheck'} if ($cgivars{'taxacheck'});
my $cutoff           = $cgivars{'cutoff'};              # Cutoff length
my $orig_stringency  = $cgivars{'stringency'};          # Stringency filter
my $bpcutoff_low     = $cgivars{'bpcutoff_low'};        # shortest acceptable fragment
my $bpcutoff_high    = $cgivars{'bpcutoff_high'};       # longest acceptable fragment
my $mismatches       = $cgivars{'mismatches'};          # max allowable group mismatches
my $matchlimit       = $cgivars{'matchlimit'};          # max matches to print, <1000
my $rdpinput         = $cgivars{'rdpinput'};            # output file from RDP classifier
( $matchlimit = 1000 ) if ( $matchlimit > 1000 );
copy("$enzymeFile","./$timestring/$enzymeFile");

my $stringency = ();
my $redo = 0;
if ($orig_stringency < 0) {
 $redo++;
 $stringency = $repeats * 0.1;
}
else {$stringency = ($orig_stringency / 100) + $repeats * 0.1;}

##### options
unless ($taxacheck) {
print INDEX "<table border = 1 cellspacing = 3 cellpadding = 2>\n";
print INDEX "<tr><td><h3>OPTIONS</h3></td><td></td></tr>\n";
print INDEX "<tr><td>Taxonomic Rank Delimiter</td><td>$splitter</td></tr>\n";
print INDEX "<tr><td>Taxonomic Rank</td><td>" . ($splitnum + 1) . "</td></tr>\n";
print INDEX "<tr><td>Cutoff</td><td>$cutoff</td></tr>\n";
print INDEX "<tr><td>Min Fragment Lengths</td><td>$bpcutoff_low</td></tr>\n";
print INDEX "<tr><td>Max Fragment Lengths</td><td>$bpcutoff_high</td></tr>\n";
print INDEX "<tr><td>Stringency</td><td>" . ($stringency * 100) . "</td></tr>\n";
print INDEX "<tr><td>Max Missing Group Combinations</td><td>$mismatches</td></tr>\n";
print INDEX "<tr><td>Max Matches Returned</td><td>$matchlimit</td></tr>\n";
# print INDEX "</table>\n";
print INDEX "<tr><td><h3>RESULTS</h3></td><td></td></tr>\n";


##### warnings
print INDEX "<tr><td><font color=red>WARNING</td><td>NO FASTA FILE DETECTED</font></td></tr>" and exit unless $fastafile;
print INDEX "<tr><td><font color=red>WARNING</td><td><em>Cutoff</em> cannot be greater than or equal to <em>Min Fragment Length</em></font></td></tr>" and exit if ( $bpcutoff_low <= $cutoff );
}
##### Which enzymes were chosen
my $seenenz = 0;
my @enzymesList = ();
$seenenz++ if length($enzymeFileCustom) != 0;

if ( length($enzyme_list) != 0 ) {
    if ( $enzyme_list =~ m/allrebase/ ) { push( @enzymesList, 'allrebase' ); $seenenz++ }
    elsif ( $enzyme_list =~ m/norebase/ ) { push( @enzymesList, 'norebase' ) }
    else { @enzymesList = split( "\0", $enzyme_list ); $seenenz++ }
}

print INDEX "<tr><td><font color=red>WARNING</td><td>NO ENZYMES CHOSEN</font></td></tr>" and exit if ( $seenenz == 0 );

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
my $fragtable      = "./$timestring/fragtable.html";
my $enztable       = "./$timestring/enztable.html";
my $matrixfig      = "./$timestring/matrixfig.png";
my $ndgfile        = ();                           # debugging file to list non-discriminated groups

##### widely used hashes
my %success = ();
my %fragments = ();
my %fragsave = ();
my %list = ();
my %partlist = ();
my %cutlength = ();
my %matchlength = ();
my %reverseEnzymes = ();
my %isoschiz = ();
my %goodcombo = ();
my %usercombo = ();
my %enzymes = ();
my %customenzymes = ();
my %userfinalsave = ();
my %finalsave = ();
my %enzymeSites = ();
my %taxasubs = ();

##### variables that are printed out of scope
my $listcount = ();
my $partcount = ();
my @nocutters = ();
my @justEnz = ();
my $combonumber     = 0;
my $usercombonumber = 0;
my $usermatchlimit  = $matchlimit;
my @taxasplit = ();
my @matrixCounter = ();
my $maxsub = 0;

##### Set up Sequences
my @renamed   = ();
@renamed = readRDPCfile( ( split( /\n/, $rdpinput ) ) );    # parse RDP-classifier output
my %sequences = ();
%sequences = readFASTAfile();                                 # parse FASTA sequences

my @subs = ();
@subs = checkGroups(%sequences);                              # make the groups


if ($taxacheck) {
print INDEX "<h2>This is what REPK thinks your taxonomic ranks are, if you expected something else please consult <a href=\"http://code.google.com/p/repk/wiki/Manual\">The Manual</a> to find out how to amend them.</h2>\n";
print INDEX "<table border = 1 cellspacing = 3 cellpadding = 2>";
print INDEX "<tr><td>Sequence Name</td>";
for (my $i; $i< $maxsub; $i++) { print INDEX "<td>taxa rank " . ($i+1) . "</td>";}
print INDEX "</tr>\n";

foreach my $seqkey ( sort { $a cmp $b } keys %sequences ) {
print INDEX "<tr><td>$seqkey</td><td>" . join("</td><td>", @{$taxasubs{$seqkey}}) . "</td></tr>\n";
}
print INDEX "</table>";
exit;
}

my $numgroups = scalar(@subs);                                   # how many groups
my $numcombos = $numgroups * ( $numgroups - 1 ) / 2;             # how many combos of those groups
my $maxfail   = $numcombos - int( $numcombos * $stringency );    # how many to pass stringency

my @grpgrp    = ();
@grpgrp = groupGroups(@subs);                              # pairwise group combinations
my @revgrpgrp = reverse @grpgrp;                                 # reverses order of list

##### Print progress on Sequences
if ( $rdpinput =~ m/Classifier/g ) {
    open( RENAMED, ">$renamedfile" ) or die "Couldn't open $renamedfile: $!";
    foreach my $seqkey ( sort { $a cmp $b } keys %sequences ) {
        print RENAMED ">$seqkey\n$sequences{$seqkey}\n";
    }
    close RENAMED or die "Couldn't close $renamedfile: $!";
    print INDEX "<tr><td><a href=\"renamed.fasta\">renamed.fasta</a></td><td>Input sequences renamed by RDP-Classifier</td></tr>\n";
}




# print INDEX "<tr><td>Pairwise group combinations</td><td>$numcombos</td></tr>\n";
print INDEX "<tr><td>Sequence Groups</td><td>$numgroups groups (" . join(', ',@subs) . ") have $numcombos possible combinations</td></tr>";
# print INDEX "<tr><td>" . (int( $numcombos * $stringency )) . "</td><td>Minimum number of group discriminations allowed by stringency</td></tr>\n";

##### Set up Enzymes
# If 'all enzymes' is chosen
if ( $enzymesList[0] eq 'allrebase' ) {
    %enzymes = readEnzfile($enzymeFile);
    print INDEX "<tr><td><a href=\"$enzymeFile\">$enzymeFile</a></td><td>List containing all " . ( scalar keys %enzymes ) . " standard enzymes</td></tr>\n";
}

# If 'all enzymes' is not chosen but others are
elsif ( $enzymesList[1] ) {
    %enzymes = readEnzchosen(@enzymesList);
    foreach (@enzymesList) { my ( $just, $junk ) = split( '\t', $_ ); push( @justEnz, $just ); }
    print INDEX "<tr><td><a href=\"custom.txt\">custom.txt</a></td><td>List containing " . scalar(@justEnz) . " user-selected enzymes</td></tr>\n";
}

if ($enzymeFileCustom) {
    %customenzymes = readEnzcustom($enzymeFileCustom);
    print INDEX "<tr><td><a href=\"custom.txt\">custom.txt</a></td><td>List containing " . ( scalar keys %customenzymes ) . " user-input enzymes</td></tr>\n";
    @enzymes{ keys %customenzymes } = values %customenzymes;    # get

}

##### print custom enzymes
open( CUSTOM, ">$customfile" ) or die "Couldn't open $customfile: $!";
print CUSTOM join( "\n", ( @justEnz, ( keys %customenzymes ) ) ) . "\n";
close CUSTOM or die "Couldn't close $customfile: $!";

#######################
####### PART 1 ########  In which the fragments are cut with the enzymes, the terminal fragment lengths are saved,
#######################  the groups are differentiated, evil is vanquished, and peace is brought to the entire galaxy.

print "\n<br><br>Stringency = " . ($stringency * 100) . " on repeat $repeats";
print "\n<br>Enzymes left to test: ";
my $enzhashsize = keys(%enzymes);
ENZYME: foreach my $RE ( sort ( keys %enzymes ) ) {

    $enzhashsize--;
    if ( ( $enzhashsize % 10 ) == 0 ) { print "$enzhashsize " }
    else { print "."}

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
        my $forwardFragment = ();
        my $reverseFragment = ();


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

#DEBUG1
            if ( $seengroup le $seenfail ) { push @{ $list{$RE} }, "$seengroup-$seenfail"; }
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
# print INDEX "PROCESSING: Enzymes successfully tested.";


##############################
##### PRINT REULTS TO FILE   #
##############################

##### print to file the enzymes that cut and distinguished
if ($ndgfile) {
    open( NDG, ">$ndgfile" ) or die "Couldn't open $ndgfile: $!";
    print NDG "#NDG\tEnzyme\tNon-discriminated groups";

    $listcount = keys %list;
    $partcount = keys %partlist;
    foreach my $okEnz ( sort ( keys %partlist ) ) {
        my $misscount = @{ $partlist{$okEnz} };
        print NDG "\n$misscount\t$okEnz\t$enzymes{$okEnz}\t" . "@{$partlist{$okEnz}}";
    }
    close NDG or die "Couldn't close $ndgfile: $!";
}

if ($yescutsfile) {
    open( YESCUTS, ">$yescutsfile" ) or die "Couldn't open $yescutsfile: $!";

    $listcount = keys %list;
    $partcount = keys %partlist;
    foreach my $okEnz ( sort ( keys %partlist ) ) {
        my $misscount = @{ $partlist{$okEnz} };
        print YESCUTS "\n$okEnz\t$enzymeSites{$okEnz}";
    }
    close YESCUTS or die "Couldn't close $yescutsfile: $!";
}


# print INDEX "<tr><td>$partcount</td><td>Enzymes that differentiated >" . ( 100 * $stringency ) . "\% (" . ( $numcombos - $maxfail ) . ") of the group combinations</td></tr>";
# print INDEX "<tr><td>" . eval { $listcount - $partcount } . "</td><td>Enzymes that did not differentiate.\n";
print INDEX "<tr><td><a href=\"yes_cuts.txt\">yes_cuts.txt</a></td><td>Enzymes that successfully passed stringency and cut sequences ($partcount)</td></tr>\n";




# ##### print fragment lengths to file
# if ($fragfile) {
#     open( FRAGS, ">$fragfile" ) or die "Couldn't open $fragfile: $!";
#     foreach my $taxaname ( sort ( keys %sequences ) ) { print FRAGS "\t$taxaname" }
#     print FRAGS "\n";
#   PRINTZYME: foreach my $printzyme ( sort ( keys %fragsave ) ) {
#         print FRAGS "$printzyme\t" . ( join "\t", @{ $fragsave{$printzyme} } ) . "\n";
#     }
#     close FRAGS or die "Couldn't close $fragfile: $!";
# }



###### print translocated fragment lengths to file
##### no good because stupid spreadsheets can only handle 256 columns!!

if ($fragfile) {

    open( FRAGS, ">$fragfile" ) or die "Couldn't open $fragfile: $!";

 foreach my $printzyme ( sort ( keys %fragsave ) ) { print FRAGS "\t$printzyme"}
    print FRAGS "\n";
my $printcount = 0;
 foreach my $taxaname ( sort ( keys %sequences ) ) {
      print FRAGS "\n$taxaname";
   foreach my $printzyme ( sort ( keys %fragsave ) ) {
      print FRAGS "\t@{$fragsave{$printzyme}}[$printcount]";
#         print FRAGS "$taxaname\t" . ( join "\t", @{ $fragsave{$printzyme} } ) . "\n";
     
    }
      $printcount++;

}
    close FRAGS or die "Couldn't close $fragfile: $!";
}

# # ##### print translocated fragment lengths to html table


if ($fragfile) {

    open( FRAGS, ">$fragtable" ) or die "Couldn't open $fragtable: $!";
print FRAGS <<EOF;
<html>
<body>
<h2>Terminal Fragment Lengths</h2>
download as tab-delimited <a href="./fragfile.csv">text file</a> for spreadsheet applications
<br>(note this file is transposed to allow for the 256-column limit in common spreadsheet programs)
<br>
<table border = "1" cellspacing = "3" cellpadding = "2">
EOF

print FRAGS "\n<tr><td></td>";
 foreach my $printzyme ( sort ( keys %fragsave ) ) { print FRAGS "\t<td>$printzyme</td>"}
print FRAGS "\n</tr>";

my $printcount = 0;
 foreach my $taxaname ( sort ( keys %sequences ) ) {
print FRAGS "\n<tr>";
      print FRAGS "\n<td>$taxaname</td>";
   foreach my $printzyme ( sort ( keys %fragsave ) ) {
      print FRAGS "\n\t<td>@{$fragsave{$printzyme}}[$printcount]</td>";
#         print FRAGS "$taxaname\t" . ( join "\t", @{ $fragsave{$printzyme} } ) . "\n";
    }
      $printcount++;
print FRAGS "\n</tr>";
}
print FRAGS "\n</table> </body> </html>";
    close FRAGS or die "Couldn't close $fragfile: $!";
}










print INDEX "<tr><td><a href=\"fragtable.html\">fragtable.html</a></td><td>Terminal restriction fragment lengths of all sequences cut with all selected enzymes</td></tr>\n";

##### print to file the loser enzymes that didn't cut or didn't distinguish at all
if ($failfile) {
    open( FAILURE, ">$failfile" ) or die "Couldn't open $failfile: $!";
    print FAILURE "The following enzymes did not cut/distinguish in the region between $bpcutoff_low and $bpcutoff_high\n" . join "\n", ( sort @nocutters );
    close FAILURE or die "Couldn't open $failfile: $!";
}
print INDEX "<tr><td><a href=\"no_cuts.txt\">no_cuts.txt</a></td><td>List of the " . scalar @nocutters . " enzymes that did not cut or distinguish any groups in the region between $bpcutoff_low and $bpcutoff_high bp.</td></tr>\n";

##### print to file the super-enzymes that singlehandedly distinguish all groups
if ($successfile) {
    open( SUCCESS, ">$successfile" ) or die "Couldn't open $successfile: $!";
    print SUCCESS "Enzyme\tcutsite\t (Success with fragments < $cutoff bp apart)";
    foreach my $successEnz ( keys %success ) {
        print SUCCESS "\n$successEnz\t$enzymeSites{$successEnz}"; #$success{$successEnz}
    }
    close SUCCESS or die "Couldn't close $successfile: $!";
}

if ( scalar keys %success ) {
print INDEX "<tr><td><a href=\"success.txt\">success.txt</a></td><td>List of the " . ( scalar keys %success ) . " enzymes that singly distinguish all the sequence groups.</td></tr>\n";
    HTMLclose();
}

#######################
####### PART 2 ########  In which the qualities of each enzyme is compared with each other enzyme and
#######################  the best group of four enzymes comes out victorious.

my %zlist = ();
%zlist = matrixMaker(%partlist);

# #######print enzmatrix
# {
# open( MATRIX, ">$matrixfile" ) or die "Couldn't open $matrixfile: $!";
# print MATRIX "\t" . ( join "\t", @revgrpgrp ) . "\n";
# 
# 
# foreach my $enzgroup ( sort ( keys %partlist ) ) {
# 
#   my (@matrixprint) = join "\t", @{ $zlist{$enzgroup} };
#   print MATRIX "$enzgroup\t@matrixprint\n";
# }
# 
# close MATRIX or die "Couldn't close $matrixfile: $!";
# }

#######print enzmatrix transposed
{
open( MATRIX, ">$matrixfile" ) or die "Couldn't open $matrixfile: $!";
 foreach my $enzgroup ( sort ( keys %partlist ) ) {
  print MATRIX "\t$enzgroup";
 }

for (my $printcounter = 0; $printcounter < @revgrpgrp; $printcounter++) {
print MATRIX "\n@revgrpgrp[$printcounter]";
 foreach my $enzgroup ( sort ( keys %partlist ) ) {
   print MATRIX "\t@{ $zlist{$enzgroup} }[$printcounter]";
 }
}
close MATRIX or die "Couldn't close $matrixfile: $!";

}


#######print enzmatrix HTML table
{
open( MATRIX, ">$enztable" ) or die "Couldn't open $enztable: $!";
print MATRIX <<EOF;
<html>
<body>
<h2>Enzyme Match Matrix</h2>
Download as tab-delimited <a href="./enzmatrix.csv">text file</a> for spreadsheet applications.
<br>*** indicates group combinations that cannot be discriminated by any selected enzymes.
<table border = "1" cellspacing = "3" cellpadding = "2">
EOF



print MATRIX "\n<tr> <td>Group Combination</td>";
 foreach my $enzgroup ( sort ( keys %partlist ) ) {
  print MATRIX "<td>$enzgroup</td>";
 }
print MATRIX "<td>Group Combination</td><td>SUM</td></td> </tr>\n";

for (my $printcounter = 0; $printcounter < @revgrpgrp; $printcounter++) {
 my $total = 0;
 print MATRIX "<tr><td>@revgrpgrp[$printcounter]</td>";
 foreach my $enzgroup ( sort ( keys %partlist ) ) {
   print MATRIX "<td>@{ $zlist{$enzgroup} }[$printcounter]</td>";
   $total += $zlist{$enzgroup}[$printcounter];
 }
print MATRIX "<td>@revgrpgrp[$printcounter]</td><td>$total</td>";

if ($total == 0 ) {print MATRIX "<td>***</td></tr>\n";}
else {print MATRIX "</tr>\n";}
}

print MATRIX <<EOF;
</table>
</body></html>
EOF

close MATRIX or die "Couldn't close $matrixfile: $!";
}

system("$perl_location $prettymatrix_location $repk_location$matrixfile > $repk_location$matrixfig");




########find out if there are any group combinations that CANNOT be distinguished

# foreach my $row (@matrixCounter) {print INDEX join ("\n",$matrixCounter[$row]) . "peace<br>"}
# no strict;
#     for my $aref ( @matrixCounter ) {
#         print INDEX "\t [ @$aref ],\n<br>";
#     }

    my @newMatrixCounter = ();
#     for (my $startx = my $x = 0; $x <= $matrixCounter[0]; $x++) {
for (my $startx = my $x = 0; $x <= 999; $x++) {
#         for (my $starty = my $y = 0; $y <= $#matrixCounter; $y++) {
for (my $starty = my $y = 0; $y <= $numcombos-1; $y++) {
            $newMatrixCounter[$y - $starty] += $matrixCounter[$x][$y];
# 	print $newMatrixCounter[$startx][$y - $starty];
        }
    }

#     for my $aref ( @newMatrixCounter ) {
#         print INDEX "[ @$aref ]";
#     }

# print INDEX "[ @newMatrixCounter]";

my $matrixFound = grep (/^0$/, @newMatrixCounter);
# print "matrixfound = $matrixFound";
if ($matrixFound > $mismatches)  {
print INDEX "<tr><td><font color=red>WARNING</font></td><td>There are $matrixFound pairwise group combinations that cannot be differentiated by any enzymes.  Please go back and amend your groups, or increase the <i>Max Missing Group Combinations</i> option to at least $matrixFound.  Examine <a href=\"enztable.html\">enztable.html</a> to identify the undifferentiated groups</td></tr>\n";
}

print INDEX "<tr><td><a href=\"enztable.html\">enztable.html</a></td><td>Enzyme match matrix</td></tr>\n";
print INDEX "<tr><td><a href=\"matrixfig.png\">matrixfig.png</td><td>Visual enzyme match matrix figure</td></tr>\n";

print "\n<br>Enzyme groups left to test ";

my %rev_groups = (); 
%rev_groups = revMatrix(%zlist);

my $countprint  = 1;
my $countprint2 = 1;
my $hashsize    = keys(%rev_groups);

my @keylist = ();
for my $k1 ( sort { $a cmp $b } keys %rev_groups ) {
    push( @keylist, $k1 );
}

KEYLOOP: for my $i ( 0 .. $#keylist - 3 ) {
    $hashsize--;
    ( ( $hashsize % 5 ) == 0 ) ? print "$hashsize " : print ".";

    for my $j ( $i + 1 .. $#keylist - 2 ) {
        for my $k ( $j + 1 .. $#keylist - 1 ) {
            for my $l ( $k + 1 .. $#keylist ) {
                my $newcompare = ( "$keylist[$i]" | "$keylist[$j]" | "$keylist[$k]" | "$keylist[$l]" );
                my $newcount   = ( $newcompare =~ tr/1// );

                if ( $newcount == $numcombos ) {
                    $combonumber++;
                    if ( $combonumber > 9999 ) { 
                         if (($stringency < 0.9) and ($redo > 0)) {
                           $repeats++;
                            goto BEGINNING;}
                         else {print INDEX "<tr><td><font color=red>WARNING</td><td>Too many results, go back and try raising the Stringency or lowering the Max Missing Groups options</font></td></tr>" and HTMLclose() }
                    }
                    push @{ $goodcombo{$combonumber} }, ( $rev_groups{ $keylist[$i] }[0], $rev_groups{ $keylist[$j] }[0], $rev_groups{ $keylist[$k] }[0], $rev_groups{ $keylist[$l] }[0] );
                    $finalsave{$combonumber} = ( ( $keylist[$i] ) =~ tr/1// ) + ( ( $keylist[$j] ) =~ tr/1// ) + ( ( $keylist[$k] ) =~ tr/1// ) + ( ( $keylist[$l] ) =~ tr/1// );
                }

                elsif ( $newcount >= ( $numcombos - $mismatches ) ) {
                    $usercombonumber++;
                    if ( $usercombonumber > 9999 ) {
                         if (($stringency < 0.9) and ($redo > 0)) {
                           $repeats++;
                            goto BEGINNING;}
                         else {print INDEX "<tr><td><font color=red>WARNING</td><td>Too many results, go back and try raising the Stringency or lowering the Max Missing Groups options</font></td></tr>" and HTMLclose() }
                    }

                    push @{ $usercombo{$usercombonumber} }, ( $rev_groups{ $keylist[$i] }[0], $rev_groups{ $keylist[$j] }[0], $rev_groups{ $keylist[$k] }[0], $rev_groups{ $keylist[$l] }[0] );
                    $userfinalsave{$usercombonumber} = ( ( $keylist[$i] ) =~ tr/1// ) + ( ( $keylist[$j] ) =~ tr/1// ) + ( ( $keylist[$k] ) =~ tr/1// ) + ( ( $keylist[$l] ) =~ tr/1// );
                }
            }
        }
    }
}

if ($outfile) {
    open( OUT, ">$outfile" ) or die "Couldn't open $outfile: $!";
    print OUT "ALL ENZYME GROUPS\nEnzymes in the same group differentiate the same sequences\nNeoschizomers are in brackets\n\nGroup Number\tGroup Members\n\n";
    my %tetragroups = ();
    my $matchnumber = 0;
    my $tetracount = 0;
    my $tally = ();
    foreach my $printkey ( keys %rev_groups ) {
    my @tetarray = ();
        $tetracount++;
        #####to print out isoschizomers
        foreach my $printkey3 ( @{ $rev_groups{$printkey} } ) {
            my $tempkey = "[$printkey3 " . join(" ",@{ $isoschiz{$printkey3} } ) . "]";

            push @tetarray, $tempkey;
            push @{ $tetragroups{$tetracount} }, ( $printkey3, @{ $isoschiz{$printkey3} } );
        }
        print OUT "$tetracount\t@tetarray\n";
    }

    ##### a mess follows, but it works
    print OUT <<EOF ;


SUCCESSFUL ENZYME SETS
Each set of enzyme groups below (see above for group members) distinguishes all of the input sequences
Score	Set Members
EOF
    ##### cycle through all combos to order them by highest score
    foreach my $printkey2 ( sort { $finalsave{$b} <=> $finalsave{$a} } ( keys %finalsave ) ) {    ##### keys are combonumbers, values are scalars
# print INDEX "<br>$printkey2 -- $finalsave{$printkey2} -- @{ $goodcombo{$printkey2}}";
        if ( $matchnumber > $matchlimit ) { last }
        $matchnumber++;

        printf OUT "\n%.2f", ( $finalsave{$printkey2} / $numcombos );

        ##### cycle through all four enzymes in goodcombo
        foreach my $printkey4 ( @{ $goodcombo{$printkey2} } ) {                ##### cycling through 4 enzymes
            foreach my $printkey5 ( keys %tetragroups ) {                                         ##### keys are tetracount numbers, values are enzymes+isoschizimers
                my @results = grep( /^$printkey4$/, @{ $tetragroups{$printkey5} } );
#                 if (@results) { print OUT "\t$printkey5 ($printkey4)"; }
		  if (@results) { print OUT "\t$printkey5"; $tally .= ("\t$printkey5\t")}
            }
        }
    }

  print OUT "\n\nQUICK OVERVIEW\nRelative Frequency of Enzyme Groups in above sets";
 foreach my $histkey ( sort {$a <=> $b} keys %tetragroups ) {
  my $totaln = $tally =~ tr/\t//;
$totaln = $totaln/2;
($totaln = 1) if $totaln ==0; #if no successful enzymes

  my $n = 0;
  while ($tally =~ /\t$histkey\t/g) { $n++ }
#   print OUT "\n$histkey: $n";
  printf OUT "\n%2.0f", ($histkey);
  my $newn = "-" x (100*$n/$totaln);
  print OUT $newn;
 }

    close OUT or die "Couldn't close $outfile: $!";
}

if ( $combonumber < 100 ) { $matchlimit = $combonumber }
if ( $combonumber > 0 ) {
    print INDEX "<tr><td><a href=\"finalout.txt\">finalout.txt</a></td><td>The top $matchlimit combinations of successful enzyme sets</td></tr>\n";
}
else { print INDEX "<tr><td><font color = red>WARNING</td><td>REPK found NO completely successful enzyme groups, please try again with different parameters (e.g. try decreasing the Stringency or raising the Max Missing Groups options).</font></td></tr>"; 


}





if ($missingoutfile) {
    open( USEROUT, ">$missingoutfile" ) or die "Couldn't open $missingoutfile: $!";
    print USEROUT "ALL ENZYME GROUPS\nEnzymes in the same group differentiate the same sequences\nNeoschizomers are in brackets\n\nGroup Number\tGroup Members\n\n";
    my %tetragroups = ();
    my $matchnumber = 0;
    my $tetracount = ();
    my $tally = ();
    foreach my $printkey ( keys %rev_groups ) {
my @tetarray;

        $tetracount++;
        #####to print out isoschizimers
        foreach my $printkey3 ( @{ $rev_groups{$printkey} } ) {
my $tempkey = "[$printkey3 " . join(" ",@{ $isoschiz{$printkey3} } ) . "]";
push @tetarray, $tempkey;
            push @{ $tetragroups{$tetracount} }, ( $printkey3, @{ $isoschiz{$printkey3} } );
        }

        print USEROUT "$tetracount\t@tetarray\n";

    }
    ##### a mess follows, but it works
    print USEROUT <<EOF ;


PARTIALLY SUCCESSFUL ENZYME SETS
Each set of enzyme groups below (see above for group members) distinguishes some of the input sequences
Score	Set Members
EOF

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
		if (@results) { print USEROUT "\t$printkey5"; $tally .= ("\t$printkey5\t")}
            }
        }
    }

  print USEROUT "\n\nQUICK OVERVIEW\nRelative Frequency of Enzyme Groups in above sets";
 foreach my $histkey ( sort {$a <=> $b} keys %tetragroups ) {
  my $totaln = $tally =~ tr/\t//;
$totaln = $totaln/2;
($totaln = 1) if $totaln ==0; #if no successful enzymes

  my $n = 0;
  while ($tally =~ /\t$histkey\t/g) { $n++ }
#   print OUT "\n$histkey: $n";
  printf USEROUT "\n%2.0f", ($histkey);
  my $newn = "-" x (100*$n/$totaln);
  print USEROUT $newn;
 }


    close USEROUT or die "Couldn't close $outfile: $!";
}

if ( $usercombonumber < 100 ) { $usermatchlimit = $usercombonumber }
if ($mismatches) {
    if ( $usercombonumber > 0 ) {
print INDEX "<tr><td><a href=\"missingout.txt\">missingout.txt</a></td><td>The top $usermatchlimit combinations of PARTIALLY successful enzyme sets, allowing for $mismatches mismatches</td></tr>\n";
    }
else { print INDEX "<tr><td><font color = red>WARNING</td><td>REPK found NO PARTIALLY successful enzyme groups, please try again with different parameters (e.g. try increasing the Max Missing Groups option).</font></td></tr>"; }
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
# 

print <<EOF ;
<html><head>
<META HTTP-EQUIV="Refresh" CONTENT="0; URL=$cgi_location$timestring/results.html">
</head><body>
<h1>PLEASE WAIT</h1>
<h2>Unique run ID and URL: <a href="./$timestring">$tempstring</a></h2>
</body></html>
EOF

}

sub printResults {


    # Finally, print out the complete HTML response page
    print INDEX <<EOF ;
<html>
<head><title>REPK Results</title></head>
<body>
<h2>REPK Results</h2>
<h3>Please see <a href="http://code.google.com/p/repk/wiki/Manual"> The Manual</a> for complete explanations of these output files</h3>
<font size=+1>
All files are accessible for at least 48 hours at <a href="$cgi_location$timestring">this URL</a>
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
    my %fileHash = (); %fileHash = readEnzymeText(@lines);
    return %fileHash;
}

sub readEnzcustom {
    my $line       = shift;
    my @lines      = split( /\n/, $line );
    my %customHash = (); %customHash = readEnzymeText(@lines);
    return %customHash;
}

sub readEnzchosen {
    my %chosenAllHash = (); %chosenAllHash = readEnzymeText(@_);
}

# subroutine readEnzymeCustom takes a CGI textarea input and returns
# a hash with enzyme name as a key and the recognition sequence as a value
# convert recognition sequence from IUPAC ambiguity code to regex

sub readEnzymeText {
    my %enzymeSub = ();
    foreach (@_) {
s/\r\n$/\n/;
s/\r$/\n/;
chomp;

        my @REarray = split(/\s+/);
        my $RE      = shift(@REarray);
        my $site    = shift(@REarray);
	chomp $site;
        shift(@REarray);    #kill the ISO
        @{ $isoschiz{$RE} } = @REarray;

	$enzymeSites{$RE} = $site;

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
    my @renarray = ();
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

    my $fastataxa = ();
    my %sequences = ();

    foreach ( split( "^", $fastafile ) ) {
        my $line = sprintf( "%s", $_ );
        $line =~ s/\r\n$/\n/;    # translate from DOS/win
        chomp $line;

        if ( $line =~ m/^\s$/ ) { next }
#DEBUG4 rem . $'
        elsif ( $line =~ m/FASTASTART/ ) { (@renamed) ? ( $fastataxa = shift(@renamed) . $splitter . $' ) : ( $fastataxa = $' )}
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
        print INDEX "<tr><td><font color=red>WARNING</td><td>Taxonomic rank too great -- lower to " . scalar(@taxasplit) . " or less</font></td></tr>" and exit if ( scalar(@taxasplit) < $splitnum );
#DEBUG < to <=
        my $grpid = $taxasplit[$splitnum];
	$maxsub = scalar(@taxasplit) if scalar(@taxasplit) > $maxsub;
	push @{ $taxasubs{$groupie} }, @taxasplit;
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
#DEBUG2 remove -
            push @grpgrp, "$sub1-$sub2";
        }
    }

    return @grpgrp;
}

##### subroutine matrixMaker makes the 1/0 enzyme cut matrix
#####

sub matrixMaker {

my $counter = 0;
    foreach my $enzgroup ( sort ( keys %partlist ) ) {

        push @{ $zlist{$enzgroup} }, ("1") x $numcombos; #fill the matrix with ones
        my @zrow = @{ $partlist{$enzgroup} };

        foreach my $zval (@zrow) {
            next unless $zval;
            my $countup   = 0;
            my $countdown = -1;
            foreach ( @{ $zlist{$enzgroup} } ) {
                next unless defined;
                $zlist{$enzgroup}[$countdown] = "0" if $zval eq $grpgrp[$countup]; #then turn to a zero if no match
                $countup++;
                $countdown--;
            }
        }

	$matrixCounter[$counter] = [ @{ $zlist{$enzgroup} } ] ;
	$counter++;

    }
    return %zlist;
}



##### subroutine revMatrix
##### has to do with increasing the efficiency by only running enzymes that have different cut patterns
#####
sub revMatrix {
    my %binary = ();
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
    my $in = ();
my %in = ();
    my $name = ();
my $value = ();

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
        if ($name !~ m/fasta/) {$value =~ s/[^-a-zA-Z0-9_\^.;\n\0\t]//g}
        elsif ($name =~ m/fasta/) {
# 		$value =~ s/^\>|\n\>/^\nFASTASTART/g;
		$value =~ s/^\>/^FASTASTART/g;
		$value =~ s/\n\>/\nFASTASTART/g;
		$value =~ s/ +/_/g;
		$value =~ s/[^-a-zA-Z0-9_\^.;\n\0\t]//g;
# # # DEBUG3 
# $value =`
                }

        $in{$name} .= "\0" if defined( $in{$name} );      # concatenate multiple vars
        $in{$name} .= $value;
    }

    return %in;

}

# Die, outputting HTML error page
# If no $title, use a default title
sub HTMLdie {
    my $msg = ();
my $title = @_;
    $title = "CGI Error" if $title eq '';
    print INDEX <<EOF ;
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

close INDEX or die "Couldn't close results.html: $!";
    exit;
}

sub HTMLclose {

    # Print close of HTML file
    print INDEX <<EOF ;
</table>
</font>
<h3>DONE</h3>
</body>
</html>
EOF

close INDEX or die "Couldn't close results.html: $!";
    exit;
}
