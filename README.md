# autodrive

## Setup Eigen
### 1. Download source code
```
git clone https://gitlab.com/libeigen/eigen.git
```
### 2. Include the directory
```
include_directories(/path/to/eigen)
```

## Setup OSQP
### 1. Download source code
```
git clone --recursive https://github.com/osqp/osqp
```
### 2. Create `build` directory
```
cd osqp
mkdir build
cd build
```
### 3. Create Makefiles and Compile OSQP and install
```
cmake -G "Unix Makefiles" ..
cmake --build .
sudo make install
```

## Setup IPOPT
[IPOPT doc](https://coin-or.github.io/Ipopt/INTERFACES.html)
### 1. Obtain necessary tools
```
sudo apt-get install gcc g++ gfortran git patch wget pkg-config liblapack-dev libmetis-dev
```
### 2. Download, build and install dependencies
https://coin-or.github.io/Ipopt/INSTALL.html#COMPILEINSTALL

### 3. Download Ipopt source code, build and install
```
git clone https://github.com/coin-or/Ipopt.git
```
