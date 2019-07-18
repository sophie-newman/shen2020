dir=$(pwd)
redshifts=$dir/redshift_list

cd ..

echo z gamma1 err1 gamma2 err2 phi_s err3 L_s err4 dtg err5 nfree redchi

OLD_IFS=$IFS
IFS=$'\n'
let counter=0
for z in $(cat $redshifts)
do
	python lf_fitter_at_z.py $z 1
	let counter=$(($counter+1))
done
