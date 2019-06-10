sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-4.7


After installation:
sudo rm -d /usr/bin/gcc
sudo ln -s /usr/bin/gcc-4.7 /usr/bin/gcc

This way your system will keep the gcc-4.8 but use by default gcc-4.7 .

If you need to go back to 4.8 :
sudo rm -d /usr/bin/gcc
sudo ln -s /usr/bin/gcc-4.8 /usr/bin/gcc