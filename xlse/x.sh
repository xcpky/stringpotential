jrun() {
    # Run xmake first if no arguments or specific flags are provided
    # if [[ $# -eq 0 || "$1" == "--onshellT" || "$1" == "--Det" || "$1" == "--Poles" ]]; then
    #     xmake || { echo "xmake failed"; return 1; }
    # fi
	xmake
	[[ $# -eq 0 ]] && return 0;
    # Run julia with provided arguments
    julia cscript.jl "$@"
    # Optionally view PNG if a known output is expected
    case "$1" in
        --onshellT) kitten icat onshellT.png ;;
        --Det) [[ "$2" == "1" ]] && kitten icat det.png || [[ "$2" == "3" ]] && kitten icat det.png ;;
        --poles) kitten icat pole.png ;;
		--onshellT_single) kitten icat onshellT_single.png ;;
		--Detsing) kitten icat detsing.png ;;
    esac
}

_jrun_completions() {
	local cur prev
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"
	prev="${COMP_WORDS[COMP_CWORD-1]}"

	local -A desc=(
		[--onshellT]="Plot onshell T amplitude"
		[--Det]="Plot determinant of I+VG"
		[--poles]="Plot poles in complex momentum plane"
		[--onshellT_single]="Plot single channel onshell-T"
		[--Detsing]="Plot single channel determinant"
		)
	
	case "$prev" in
		--Det)
			COMPREPLY=( $(compgen -W "0 1 2 3" -- "$cur") )
			return 0
		;;
	esac

	local opts="${!desc[*]}"
	local matches=($(compgen -W "${opts}" -- "$cur"))

	if [[ ${#matches[@]} -gt 0 ]]; then
		printf "\n"
		for opt in "${matches[@]}"; do
			printf "  %-20s %s\n" "$opt" "${desc[$opt]}"
		done | column -t
		printf "\n"
		COMPREPLY=("${matches[@]}")
	fi
	

}
complete -F _jrun_completions jrun
