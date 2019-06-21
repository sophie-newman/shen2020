dir=$(pwd)
redshifts=$dir/redshift_list

cd ..

OLD_IFS=$IFS
IFS=$'\n'
let counter=0
for z in $(cat $redshifts)
do
	#echo $z
	python lf_fitter_at_z.py $z 0
	let counter=$(($counter+1))
done
