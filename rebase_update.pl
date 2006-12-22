#!/usr/bin/perl -w
use strict;

# use LWP::Simple;
# use Data::Dumper;

#my $url = 'ftp://ftp.neb.com/rebase/link_itype2';
#my $rebase = get($url);
#print "blah $rebase";
`rm link_itype2`;
`wget --user-agent="Mozilla" http://rebase.neb.com/rebase/link_itype2`;

my %Enzymes;

open(FILEIN, "< link_itype2");

while (<FILEIN>) {
 #if the line has a tab...
 if (m/\t/) {
  #then split the line by tabs
  (my $enzyme, my $prototype, my $recsite, my $junk, my $comsource) = split("\t");
  chomp ($enzyme, $prototype, $recsite, $comsource);
  if ($comsource) {
#   $recsite =~ s/\^/_/g;

  $Enzymes{$recsite}{'comsource'} .= "$comsource,"; #append the commercial sources

my $tempproto = $Enzymes{$recsite}{'prototype'} if (defined $Enzymes{$recsite}{'prototype'});

  if (length($prototype) > 0) {
   if (defined $tempproto) {
    if ($tempproto !~ m/$prototype/g) {
#      print "$tempproto ($Enzymes{$recsite}{'prototype'}) ne $prototype -- $recsite -- $enzyme\n";
     $Enzymes{$recsite}{'prototype'} .= "$prototype";
#      print "$tempproto ($Enzymes{$recsite}{'prototype'}) ne $prototype -- $recsite -- $enzyme\n";
     }
    else {}
   }
   else {
    $Enzymes{$recsite}{'prototype'} = "$prototype"; 
    $Enzymes{$recsite}{'isos'} = "$enzyme";
   }
   if (defined($Enzymes{$recsite}{'isos'})) {
    if ($Enzymes{$recsite}{'isos'} !~ m/$enzyme[\t]+/) {$Enzymes{$recsite}{'isos'} .= "\t$enzyme"}
    else {}
   }
   else {$Enzymes{$recsite}{'isos'} = "$enzyme"}
  }
  else {
   if (defined $Enzymes{$recsite}{'prototype'}) {
    if ($Enzymes{$recsite}{'prototype'} !~ m/$enzyme/g) {$Enzymes{$recsite}{'prototype'} .= "$enzyme"}
    else {}
   }
   else {$Enzymes{$recsite}{'prototype'} = "$enzyme"}
 }
}
}
}
close FILEIN or die "complete failure";

my %NewEnzymes;
foreach my $enzkey (keys %Enzymes) {
my $temptype = $Enzymes{$enzkey}{'prototype'};
$NewEnzymes{$temptype}{'isos'} = $Enzymes{$enzkey}{'isos'};
$NewEnzymes{$temptype}{'recsite'} = $enzkey;

 if (defined $Enzymes{$enzkey}{'comsource'}) {
  my @Tempcom;
  s//$Enzymes{$enzkey}{'comsource'}/;
  while (/\S/) {push(@Tempcom,$&); s/[,]*$&[,]*//g}

  $NewEnzymes{$temptype}{'comsource'} = join ('',sort {$a cmp $b} @Tempcom);
 }
 else {$NewEnzymes{$temptype}{'comsource'} = 'none'; print 'WAIT A SECOND'}

}

my %WeirdEnzymes;
foreach my $enzkey (keys %NewEnzymes) {
 if ($NewEnzymes{$enzkey}{'recsite'} =~ m/\(/) {
  $WeirdEnzymes{$enzkey} = $NewEnzymes{$enzkey};
  delete $NewEnzymes{$enzkey};
 }
}

# print Dumper %NewEnzymes;
# 
open(OUT, "> enzymes_type2.txt");
foreach my $enzkey (sort {$a cmp $b} keys %NewEnzymes) {
print OUT "$enzkey\t$NewEnzymes{$enzkey}{'recsite'}\tISO\t$NewEnzymes{$enzkey}{'isos'}\n";
}
close OUT or die "complete failure";

open(OUT, "> enzymes_weird.txt");
foreach my $enzkey (sort {$a cmp $b} keys %WeirdEnzymes) {
print OUT "$enzkey\t$WeirdEnzymes{$enzkey}{'recsite'}\tISO\t$WeirdEnzymes{$enzkey}{'isos'}\n";
}
close OUT or die "complete failure";


