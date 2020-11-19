# cp the header 
rm -r merged
mkdir merged
cd chr11
for i in `ls` 
do
	echo $i
	head -1 ${i} > ../merged/${i} 
done
cd ..


# cp the block
for chr in {1..22}
do 
	cd chr${chr}

	for i in `ls` 
	do
		echo $i
		tail -n +2 -q ${i} >> ../merged/${i}
	done
	cd ..
done
