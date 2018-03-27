mkdir -p logs

if [[ $# -le 4 ]]; then
        echo "USAGE: $0 <BED> <GENOME> <NAME> <SIZE> <HIST> <TAGDIR> [<TAGDIR> ...]"
        exit 1
fi

BED=$1
GENOME=$2
NAME=$3
SIZE=$4
HIST=$5
shift 5
annotatePeaks.pl $BED $GENOME -size ${SIZE} -hist ${HIST} -d "$@" > ${NAME}.txt 2> logs/${NAME}.err
