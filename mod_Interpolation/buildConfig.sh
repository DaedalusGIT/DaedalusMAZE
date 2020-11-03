#1--> orbit file
#2-->run file
#3--> run path
#4--> number of model files
cfgFile=common.h

cp Template/common.h . 

if [ $# -ne 4 ]
then
  echo "Usage: ./buildConfig.sh <orbitfile> <sample model file> <Run path>
  <number of model files>"
  exit 1
fi

#Get time dimensionGet time dimension
runpath=$3
runfile=$2
orbitfile=$1
numFiles=$4

timedim=$(ncdump -h $orbitfile |grep "time =" | grep -o '[0-9]\+')
#orbitPath=$(realpath --relative-to=$runpath $orbitfile)
orbitPath=$(realpath  $orbitfile)
ntimes=$(ncdump -h $runfile |head -n 5 | grep "time =" | grep -o '[0-9]\+')



#Time to SED
sed -i 's/const int DIM =;/const int DIM = '$timedim';/g' $cfgFile

sed -i "s|OrbitName =;|OrbitName ='"$orbitPath"';|g" $cfgFile

sed -i "s|NTIMES =;|NTIMES ="$ntimes";|g" $cfgFile

sed -i 's|const int BATCH =;|const int BATCH = '$timedim'/'$numFiles';|g' $cfgFile

sed -i "s/'/\"/g"  $cfgFile 

make preload
cp InterpBin $runpath 
echo "Orbit time dimension is " $timedim
echo "Relative Orbit Path is " $orbitPath
echo "Model time dimension is " $ntimes
echo "Runpath" $runpath

