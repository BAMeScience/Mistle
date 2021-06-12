
cmd="./mistle-build -n 8 -t 1 -o ./benchmarking/ -i /home/ynowatzk/data/IGC/msp/"
echo $cmd | tee -a ./mistle.log
$cmd | tee -a ./mistle.log

cmd="./mistle-search -p 10.0 -b 0.1 -t 1 -i ./benchmarking/ -s /home/ynowatzk/data/IGC/mgf/F07_label.mgf"
echo $cmd | tee -a ./mistle.log
$cmd | tee -a ./mistle.log
