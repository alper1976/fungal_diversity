SETTING VARIABLES
VSEARCH=$(which vsearch)
SPLITS=40

mkdir ITS2
mkdir 5_8S

TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)

for FILE in `ls *.fas` ; do

      
        faSplit sequence "${FILE}" $SPLITS splitted

        for f in `ls splitted*` ; do
                ITSx -i $f --save_regions 5.8S,ITS2 --complement F -t all --preserve T --partial 50 -o $f.out &

        done

        wait

        cat splitted*ITS2.full_and_partial.fasta >> "${TMP_FASTA1}"
	cat splitted*5_8S.full_and_partial.fasta >> "${TMP_FASTA2}"

        rm splitted*


# Sorting
"${VSEARCH}" --sortbysize "${TMP_FASTA1}" \
             --output ITS2/"${FILE}" > /dev/null

# Sorting
"${VSEARCH}" --sortbysize "${TMP_FASTA2}" \
             --output 5_8S/"${FILE}" > /dev/null



rm "${TMP_FASTA1}"
rm "${TMP_FASTA2}"

done 