# converts tiffs to gray scale
#!/usr/bin/bash 


for FILE in *.tif; do 
	f=${FILE##*/}
	newf="gray/${f%.tif}_gray.tif"
	echo "Processing $f"
	convert $FILE -colorspace Gray $newf
done
