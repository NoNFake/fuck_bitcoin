# Fuck bitcoin by NoNFake


## Bulid

```
sudo pacman -S tbb
```

```bash
g++ -O3 -march=native -funroll-loops -ffast-math test.cpp -o pohlig_hellman -lgmpxx -lgmp -lpthread -ltbb
```

### RUN
```
./pohlig_hellman
```



