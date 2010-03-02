<html>
<body>
<?php
 $repk_location="/Library/Webserver/Documents/cgi/repk/";
 $repk_filename="repk.cgi";
 $cgi_location="/cgi/repk/";
?>

<h2>Restriction Enzyme Picker Online, v.1.3</h2>
by Gabrielle Rocap and Eric Collins (rec3141.at.gmail.com)
<br><br>
REPK (Restriction Endonuclease Picker) finds sets of 4 commercially available restriction endonucleases 
which together uniquely differentiate designated sequence groups from a supplied FASTA format 
sequence file for use in T-RFLP.  If you have a database of known sequences from an environment, this 
program could be useful to pick some enzymes that uniquely discriminate the different groups in your database.
<br><br>
<big><strong><a href="http://code.google.com/p/repk/wiki/Manual">Read the Manual</a></big></strong>
<br>
<strong><font color="red">REPK may take HOURS when a lot of groups are used. See note below.</font></strong>
<br><br>
 <div id="formDiv">
<form action="<?php print $cgi_location.$repk_filename;?>" method="post" name="dataform" 
method="post">

<p>Enter your FASTA aligned sequences here:
<br>(Here is an example input file: <a href="<?php print $cgi_location."alignment5.txt";?>">alignment5.txt</a>)<br>
<TEXTAREA NAME="fasta" COLS=90 ROWS=10></TEXTAREA>
<br>
<p>Select the enzymes you would like to use: (<a href="<?php print $cgi_location."enzymes_type2p.txt";?>">REBASE Version 

<?php
system("cat $repk_location" . "rebase_version.txt",$version);
substr($version, 0, -1);
?>
</a>)

<br>(Hold down the [Control] key to deselect, or to select multiple enzymes)
<br>At least two enzymes must be used.
<br>NEB list of <a href="http://www.neb.com/nebecomm/tech_reference/restriction_enzymes/isoschizomers.asp" target="_new">isoschizomers</a>, the first (alphabetically) is shown below.
<select name="enzyme_list" multiple>
<option value='allrebase' SELECTED>All Commercially Available Type IIP Enzymes in REBASE</option>
<option value='norebase'>Use my Custom enzymes only</option>

<?php
$enzFile = "$repk_location" . "enzymes_type2p.txt";
$fh = fopen($enzFile, 'r');
if ($fh) {
 while (!feof($fh)) {
  $enzline = fgets($fh, 4096);
  //echo "$enzline<br>";
  list($enzsite, $isos) = split("ISO\t", $enzline);
  list($enzyme, $recsite) = split("\t", $enzsite);
  $enzyme = rtrim($enzyme);
  $newisos = rtrim(join(' ',split("\t",$isos)));
  print "<option value='$enzyme\t$recsite'>$enzyme [$recsite] [$newisos]</option>\n";

 }
}
fclose($fh);
?>

</select>
<br><a href="<?php print $cgi_location."enzymes_type2a.txt";?>">Type IIA enzymes (forward and reverse [+])</a>
<br><br>
Enter any custom enzymes you would like to use in the following format (where ^ is the cut site): <br><em>NoCut [tab] zzz^zzz</em> 
<br>
<textarea name="enzymeFileCustom" cols="20" rows="2" wrap="off"></textarea>


<!--
<p><em>Group Subset Delimiter</em>
<select name="splitter">
<option value="_" SELECTED>_ (underscore)</option>
<option value="|">| (pipe)</option>
</select>: the character used to separate groups in your FASTA sequence names
-->

<p><em>Taxonomic rank</em>: <input type="text" name="splitnum" size=4 maxlength=4 value="1">: which taxonomic group do you want separated out?
<br><input type="checkbox" name="taxacheck"><i>CHECK TAXONOMIC RANKS INSTEAD OF RUNNING REPK</i>



<p><em>Cutoff</em>: <input type="text" name="cutoff" size=4 maxlength=4 value="5">: furthest apart two fragments can be in length and still be considered the same fragments
<p><em>Min Fragment Lengths</em>: <input type="text" name="bpcutoff_low" size=4 maxlength=4 value="75">: shortest acceptable fragment (must be greater than <em>Cutoff</em>)
<p><em>Max Fragment Lengths</em>: <input type="text" name="bpcutoff_high" size=4 maxlength=4 value="900">: longest acceptable fragment
<p><em>Stringency</em>: 
<select name='stringency'>
<option value='-1' SELECTED>automatic</option>
<option value='0'>0</option>
<option value='10'>10</option>
<option value='20'>20</option>
<option value='30'>30</option>
<option value='40'>40</option>
<option value='50'>50</option>
<option value='60'>60</option>
<option value='70'>70</option>
<option value='80'>80</option>
<option value='90'>90</option>
<option value='100'>100</option>
</select>: an enzyme must distinguish MORE than this percent of groups to be acceptable
<p><em>Max Missing Group Combinations</em>: <input type="text" name="mismatches" size=4 maxlength=2 value="0">: the number of missing groups allowed
<p><em>Max Matches Returned</em>: <input type="text" name="matchlimit" size=4 maxlength=4 value="100">: the maximum number of matches to print out (max allowed 1000)

<p><input type="submit" name="form_done" value="Submit Form">
<input type="reset" value="Reset Form"><br><br>

<p><em><b>OPTIONAL</b></em> -- Enter your RDP Classifier output file here -- <em><b>OPTIONAL</b></em>
<br>(Here is an example input file: <a href="<?php print $cgi_location."rdpdownload5.txt";?>">rdpdownload5.txt</a> or generate your own here: <a href="http://rdp.cme.msu.edu/classifier/">RDP-Classifier</a>)<br>

<TEXTAREA NAME="rdpinput" COLS=90 ROWS=10></TEXTAREA>
<br>

</form>
</div>
<br><br>
<font color="red"><strong>NOTE: Some recent runs with more than 75 groups and 179 enzymes took 4 hours to complete.  Don't close the page until it says "DONE". If you close the page then the results.html file will not update but the program will continue running. As of 26/01/2010 You can check for these files by browsing to the base directory of the REPK run, e.g. http://rocaplab.ocean.washington.edu/cgi/repk/RESULTS/REPK_2010026-1323-32027e46281c7b3c/ where REPK_2010026-1323-32027e46281c7b3c is the unique run ID.
</strong></font>


<br><br>
<big><strong>CITE REPK</big></strong><br>
If you find REPK useful in your research please cite:</p><p>Collins, R. E. and G. Rocap.  2007.  REPK: an analytical web server to select restriction endonucleases for terminal restriction fragment length polymorphism analysis.  <em>Nucleic Acids Res.</em>  35 (Database issue): W58-W62; doi:10.1093/nar/gkm384   [<a href="http://nar.oxfordjournals.org/cgi/content/full/35/suppl_2/W58?ijkey=8OPX1p1IyGKNWuj&amp;keytype=ref">Free full text</a>]<br /><br /> Source code, help, and license information is available via Google Code:<br /> <a href="http://code.google.com/p/repk/">http://code.google.com/p/repk</a> <br /> <br /> Please report bugs or other problems here:<br /> <a href="http://code.google.com/p/repk/issues/list">http://code.google.com/p/repk/issues/list</a>  <br /><br />Version Changes <br>1.3 -- Added directory listing <br />1.2 -- Added Visual Enzyme Match Matrix <br />1.1 -- Added Taxonomic Rank Checker <br />1.0 -- Initial Release</p><p>&nbsp;</p><p><img src="/system/files/nsf4c.jpg" alt="NSF logo" title="NSF logo" width="100" height="100" align="left" /> </p><p>Support for the development of REPK provided by the National Science Foundation under Grant No. OPP-0327244 to J.W. Deming and OCE-0822026 and OCE-0352190 to G. Rocap, as well as a Washington<sup> </sup>Sea Grant award to J.W. Deming.</p><p> Any opinions, findings and conclusions or recomendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation (NSF) </p>
</body>
</html>