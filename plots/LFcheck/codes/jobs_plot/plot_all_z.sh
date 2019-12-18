dir=$(pwd)
redshifts=$dir/redshift_list_all

cd ..

OLD_IFS=$IFS
IFS=$'\n'
let counter=0
for z in $(cat $redshifts)
do
	python plot_allbands.py $z
	let counter=$(($counter+1))
done
