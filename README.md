# SummeStudent2019Tag-Probe
Instruction on how to run the code:

1) copy the file trackingNTuple.root
2) open root: `root -l`
3) run the code:

```
root [0] .L myReader.C++g
root [1] myReader t
root [2] t.Loop()
```