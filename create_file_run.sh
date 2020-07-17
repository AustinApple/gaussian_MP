for i in 14 21 45 63 70 72 77 81 95 
do
mkdir $i 
cd $i 
mkdir neu pos neg
cp ../${i}_neu.com neu/
cp ../${i}_pos.com pos/
cp ../${i}_neg.com neg/

cd neu 
g03shell_MP ${i}_neu
qsub g03sub.sh
cd ../ 

cd pos
g03shell_MP ${i}_pos
qsub g03sub.sh
cd ../

cd neg
g03shell_MP ${i}_neg
qsub g03sub.sh
cd ../

cd ../ 
done
