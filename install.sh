

sudo cat nice.txt


echo
echo "Checking for the dependencies . . ."

echo -n "Checking for Bowtie2: "

which bowtie2 &> /dev/null

if [ $? -eq 1 ]; then
    echo "Bowtie2 is missing! Please install Bowtie2."
    exit
else
    echo "OK"
fi


echo -n "Checking for Networkx python library: "

python2.7 -c "import networkx" &> /dev/null

if [ $? -eq 1 ]; then
    echo "Networkx is missing! Please install Networkx."
    exit
else
    echo "OK"
fi



echo -n "Checking for Numpy python library: "

python2.7 -c "import numpy" &> /dev/null

if [ $? -eq 1 ]; then
    echo "Numpy is missing! Please install Numpy."
    exit
else
    echo "OK"
fi






sudo cp -r ../ScaffMatch /usr/bin
sudo ln -fs /usr/bin/ScaffMatch/scaffmatch /usr/bin/scaffmatch

echo "Congratulations! You have successfully installed Scaffmatch on your system."
