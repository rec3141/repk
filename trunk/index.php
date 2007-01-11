<?php include("header_doctype.txt");?>
<html>
<head>

<title>REpicker</title>
</head>

<body>

<h2>Restriction Enzyme Picker Online</h2>
by Gabrielle Rocap and Eric Collins (rec3141.at.gmail.com)
<br><br>
REPK (Restriction Endonuclease Picker) finds sets of 4 commercially available restriction endonucleases 
which together uniquely differentiate designated sequence groups from a supplied FASTA format 
sequence file for use in T-RFLP.  If you have a database of known sequences from an environment, this 
program could be useful to pick some enzymes that uniquely discriminate the different groups in your database.
<br><br>
<big><strong><a href="http://code.google.com/p/repk/wiki/Manual">Read the Manual</a></big></strong>
<br><br>
 <div id="formDiv">
<form action="./repk0.2.8.cgi" method="post" name="dataform" method="post">
<p>Enter your FASTA aligned sequences here:
<br>(Here is an example input file: <a href="./alignment.4.fas">alignment.4.fas</a>)<br>
<TEXTAREA NAME="fasta" COLS=90 ROWS=10></TEXTAREA>
<br>
<p>Select the enzymes you would like to use: (<a href="./enzymes_type2.txt">see plaintext list</a>)
<br>(Hold down the [Control] key to deselect, or to select multiple enzymes)<br>
<select name="enzyme_list" multiple>
<option value='allrebase' SELECTED>All Commercially Available Enzymes in REBASE</option>
<option value='norebase'>Use my Custom enzymes only</option>

<?php
$enzFile = "./enzymes_type2.txt";
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
<br><br>
Enter any custom enzymes you would like to use in the following format (where ^ is the cut site): <br><em>NoCut zzz^zzz</em> 
<br>
<textarea name="enzymeFileCustom" cols="20" rows="2" wrap="off"></textarea>


<p><em>Group Subset Delimiter</em>
<select name="splitter">
<option value="_" SELECTED>_ (underscore)</option>
<option value="|">| (pipe)</option>
</select>: the character used to separate groups in your FASTA sequence names

<p><em>Group Subset</em>: <input type="text" name="splitnum" size=4 maxlength=4 value="1">: which group subset do you want separated out?
<p><em>Cutoff</em>: <input type="text" name="cutoff" size=4 maxlength=4 value="5">: furthest apart two fragments can be in length and still be considered the same fragments
<p><em>Min Fragment Lengths</em>: <input type="text" name="bpcutoff_low" size=4 maxlength=4 value="75">: shortest acceptable fragment (must be greater than <em>Cutoff</em>)
<p><em>Max Fragment Lengths</em>: <input type="text" name="bpcutoff_high" size=4 maxlength=4 value="900">: longest acceptable fragment
<p><em>Stringency</em>: <input type="text" name="stringency" size=4 maxlength=4 value="0">: an enzyme must distinguish MORE than this percent of groups to be acceptable
<p><em>Max Missing Groups</em>: <input type="text" name="mismatches" size=4 maxlength=2 value="0">: the number of missing groups allowed
<p><em>Max Matches Returned</em>: <input type="text" name="matchlimit" size=4 maxlength=4 value="100">: the most number of matches to print out (maximum 1000)

<p><input type="submit" name="form_done" value="Submit Form">
<br>
<p><em><b>OPTIONAL</b></em> -- Enter your RDP Classifier output file here -- <em><b>OPTIONAL</b></em>
<br>(Here is an example input file: <a href="./rdpdownload.4.txt">rdpdownload.4.txt</a> or generate your own here: <a href="http://rdp.cme.msu.edu/classifier/">RDP-Classifier</a>)<br>

<TEXTAREA NAME="rdpinput" COLS=90 ROWS=10></TEXTAREA>
<br>

</form>
</div>
<br><br>
Source code, help, and license information is available via Google Code:<br>
<a href="http://code.google.com/p/repk/">http://code.google.com/p/repk</a>
<br>
<br>
Please report bugs or other problems here:<br>
<a href="http://code.google.com/p/repk/issues/list">http://code.google.com/p/repk/issues/list</a>




</html>
</body>
