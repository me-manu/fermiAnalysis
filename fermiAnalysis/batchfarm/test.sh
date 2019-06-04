br=`bjobs -r | grep RUN | wc -l`
bp=`bjobs -p | grep PEN | wc -l`
echo $br
echo $bp
echo $(($bp + $br))

sleep=15.
i=0

while [ $i -lt 10 ]
    do
        echo $i
        while [ $(($bp + $br)) -gt 40 ]
            do

                echo "RUN $br PEND $bp"
                echo "Sleeping for $sleep seconds"
                sleep $sleep
                br=`bjobs -r | grep RUN | wc -l`
                bp=`bjobs -p | grep PEN | wc -l`
            done
        i=$(($i + 1))
    done
