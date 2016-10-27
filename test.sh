#!/bin/zsh

# scons unsets this
export LC_ALL=C

function normalize() {
    next=0
    current=( )
    step=0
    result="ERROR"
    while read line; do
        if [[ next -eq 1 ]]; then
            if [[ step -gt 0 ]]; then
                print "Step: $step"
                current=(${(@on)current})
                for answer in "${current[@]}"; do
                    print $answer
                done
                current=( )
            fi
            step=$[step+1]
            next=0
        fi
        if [[ next -eq 2 ]]; then
            answer=(${(s: :)line})
            answer=(${(@on)answer})
            current+=("${answer[*]}")
            next=0
        elif [[ "$line" =~ "^Solving..." ]]; then
            next=1
        elif [[ "$line" =~ "^Answer: " ]]; then
            next=2
        elif [[ "$line" =~ "^SATISFIABLE" ]]; then
            result="SAT"
            next=1
        elif [[ "$line" =~ "^UNSATISFIABLE" ]]; then
            result="UNSAT"
            next=1
        elif [[ "$line" =~ "^UNKNOWN" ]]; then
            result="UNKNOWN"
            next=1
        elif [[ "$line" =~ "^OPTIMUM FOUND" ]]; then
            result="OPTIMUM FOUND"
            next=1
        fi
    done
    print "$result"
}

function usage() {
    cat << EOF
Usage:
  test.sh {-h,--help,help}
  test.sh [PATH-TO-GRINGO] [PATH-TO-FOUNDED] [PATH-TO-CLINGCON] [-- CLINGO-OPTIONS]
  test.sh normalize FILE [PATH-TO-GRINGO] [PATH-TO-FOUNDED] [PATH-TO-CLINGCON] [-- CLINGO-OPTIONS]

The first invocation prints this help, the second runs all tests, and the third
takes a logic program, runs grounder, translator, and solver, and normalizes
its output.
EOF
}

if [[ $# > 0 && ( "$1" == "--help" || $1 == "-h" || $1 == "help" ) ]]; then
    usage
    exit 0
fi
wd=$(cd "$(dirname "$0")"; pwd)
norm=0
if [[ $# > 0 && "${1}" == "normalize" ]]; then
    norm=1
    shift
    if [[ $# == 0 ]]; then 
        usage
        exit 1
    fi
    file="$1"
    shift
fi
gringo="gringo"
founded="./lc2casp"
clingcon="clingcon"
if [[ $# > 0 && "$1" != "--" ]]; then
    gringo="$1"
    shift
    if [[ $# > 0 && "$1" != "--" ]]; then
        founded="$1"
        shift
        if [[ $# > 0 && "$1" != "--" ]]; then
            clingcon="$1"
            shift
        fi
    fi
fi
if [[ $# > 0 && "$1" != "--" ]] then
    usage
    exit 1
fi
[[ $# > 0 && "$1" == "--" ]] && shift
if [[ $norm == 1 ]]; then
    name=${file%.lp}
    opts=( )
    if [[ -e "$name.cmd" ]]; then
        while read line; do
            opts+=("$line")
        done < <(cat "$name.cmd")
    fi
    co=(${(s: :)opts[1]})
    fo=(${(s: :)opts[2]})
    go=(${(s: :)opts[3]})
    $gringo "$file" "${go[@]}" | $founded "${fo[@]}" | $clingcon 0 "${co[@]}" "$@" | normalize
    exit 0
else
    run=0
    fail=0
    failures=()
    for x in $wd/test/**/*.lp; do
        run=$[run+1]
        name=${x%.lp}
        opts=( )
        if [[ -e "$name.cmd" ]]; then
            while read line; do
                opts+=("$line")
            done < <(cat "$name.cmd")
        fi
        co=(${(s: :)opts[1]})
        fo=(${(s: :)opts[2]})
        go=(${(s: :)opts[3]})
        if $gringo "$x" "${go[@]}" | $founded "${fo[@]}" | $clingcon 100 "${co[@]}" "$@" | normalize | diff - "$name.sol"; then
            print -n "."
        else
            print -n "F"
            fail=$[fail+1]
            failures+=($x)
            #TODO: record failure and report at the end
        fi
    done
    print
    print
    print -n "OK ($[run-fail]/${run})"
    print
    print
    if [[ fail -gt 0 ]]; then
        print "The following tests failed:"
        for x in "${failures[@]}"; do
            print "  $x"
        done
    fi
    [[ fail -eq 0 ]]
fi
