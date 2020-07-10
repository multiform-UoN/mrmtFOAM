cat log | grep breakthrough > breakthrough.dat
echo "$(tail -n +2 breakthrough.dat)" > breakthrough.dat
python3 plot.py
