mkdir tmpzip
mv  ./document/allSNV.xlsx ./tmpzip/
mv  ./vcf/sample.all.final.vcf.gz ./tmpzip/
mv  ./vcf/sample.all.final.vcf.gz.tbi ./tmpzip/
cd tmpzip
zip $1_OriginalResult.zip  *
mv allSNV.xlsx ../document/
mv sample.all.final.vcf.gz ../vcf/
mv sample.all.final.vcf.gz.tbi ../vcf/
mv * ../
cd ..
rm -r tmpzip
