#!/bin/bash

wget http://pub.ist.ac.at/~vnk/software/QPBO-v1.3.src.tar.gz
tar -xzf QPBO-v1.3.src.tar.gz

wget http://www.f.waseda.jp/hfs/HOCR1.02.zip
unzip HOCR1.02.zip -d tmp
mv tmp/Image.h HOCR/Image.h
mv tmp/HOCR/HOCR.h HOCR/HOCR.h
mv tmp/HOCR/HOCR0.h HOCR/HOCR0.h

wget ftp://ftp.cs.man.ac.uk/pub/toby/gpc/gpct101.zip
unzip gpct101.zip -d tmp
mv tmp/gpct101/gpc.c thirdparty/gpc.c
mv tmp/gpct101/gpc.h thirdparty/gpc.h

rm QPBO-v1.3.src.tar.gz
rm HOCR1.02.zip
rm -r tmp

