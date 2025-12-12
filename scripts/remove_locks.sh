#!/usr/bin/env bash
TARGET="./data/signals/"
for i in `ls ${TARGET}|grep -P ".lock$"`;do 
    NO_LOCK=`printf ${i}|sed -e 's/.lock//'`
    echo "renaming ${i} -> ${NO_LOCK}";
    mv ${TARGET}/${i} ${TARGET}/${NO_LOCK}
done
rm ${TARGET}/STOP 