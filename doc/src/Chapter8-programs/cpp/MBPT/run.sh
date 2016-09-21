method=2

for function in 0 1 2 3 4; do
    ./MBPT 0.08 28 0 2 2 ${function}
    echo
done
