# ChIP-seq prep!

TF_NAMES=my_tf_names.txt
# note: these files 
TF_FILES=k562_tfbs_noIfng.Mnase.Ifna.txt

## -- Now for ChIP-seq! -- #
# #bigBedToBed raw_chip bedchip...
# # Need to look for ChIP-seq peaks in the DNase peaks! the peak list is $PEAK_LIST.totals ... alternately $PEAK_SEQ_WITHSIGNAL.gz should contain this information (these should agree on peaks... need to double check this! maintain consistency!)

#For obtaining the data... let's assume I have a list of all the files on the ENCODE TFBS download page, then we just basically have
cat webpage.txt | grep -i 'k562' | grep -i 'tfbs' | sed '/Ifng/d' > k562_tfbs_noIfng.txt
# this is pretty horrifying, why am i not just using python WHAT PURPOSE DOES THIS EVEN SERVE
cut -d '/' -f 15 k562_tfbs_noIfng.txt | cut -d '.' -f 3 | cut -c 9- | cut -d '_' -f 1 | cut -d '2' -f 2- | sed 's/Rep0//' > almost_Tfs.txt
# nate gave me a list of ENCODE TFs, so I can scoop those out from the master list (or the aesthetically more pleasing result of the previous command)
for i in {1..127}
do
    tf_name=`sed ''$i'q;d' TF_names.txt`
    echo $tf_name `grep -i "$tf_name" almost_TFs.txt` >> my_tfs.txt
done
# some of these aren't found and simply produce boring columns
awk 'NF-1' my_tfs.txt | uniq > to_get.txt
# the result of this seems to be 57 TFs... ... save these! ... BAD SCRIPT!
for i in {1..57}
do
        tf_ident=`sed ''$i'q;d' to_get.txt | cut -d ' ' -f 2`
             remote_file=`grep $tf_ident webpage.txt`
             tf_name=`sed ''$i'q;d' to_get.txt | cut -d ' ' -f 1`
             wget $remote_file -o $tf_name.bb
         done
         zcat $CHIP_PEAKS | bedmap --echo --indicator --delim '\t' $PEAK_LIST.totals - > $FACTOR_BOUND_PEAKS
