#!/bin/bash
TARGET="./data/signals/"
for i in `ls ${TARGET}`;do 
    echo "touching ${i}";
    touch ${TARGET}/${i}
done

