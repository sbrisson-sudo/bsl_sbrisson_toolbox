
if [ $# -lt 1 ]; then
    echo "Usage : clean_sandwich_cSEM <stock_freq|generate_bd|hsemm>"
    exit
fi



case $1 in 

    "stock_freq")
        echo "not implemented yet."
        ;;

    "generate_bd")
        rm -rfv mailleur* outfile* check_modele* fort.27 RSknfile macromesh.dat.log
        rm -rfv DATA/spem DATA/interp DATA/flt
        ;;
    
    "hsemm")
        rm -v -rf solveur* spectre.gnu time_control.log *.gnu source.gnu stop_en_run outmain* out_main Check*  *green* macromesh.dat.log fort.*
        ;;

    *)
        echo "$1 : invalid entry."
        ;;

esac

