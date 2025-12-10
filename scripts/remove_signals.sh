#!/usr/bin/env bash
TARGET="./data/signals/"
for i in `ls ${TARGET}|grep -P ".lock$"`;do 
    echo "removing ${i}";
    rm ${TARGET}/${i}
done
for i in `ls ${TARGET}|grep -P ".json$"`;do 
    echo "removing ${i}";
    rm ${TARGET}/${i}
done
rm ${TARGET}/STOP