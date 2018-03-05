rm -rf log.txt
./serial_dummy -n $1 -o log.txt
./autocorrect -s log.txt


rm -rf log.txt
./serial -n $1 -o log.txt
echo ""
./autocorrect -s log.txt
