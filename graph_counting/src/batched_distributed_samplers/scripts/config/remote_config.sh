for i in {1..31}
do
	scp ./config.sh node$i:/users/zz_y/
	ssh node$i 'cd /users/zz_y ; bash ./config.sh'
done
