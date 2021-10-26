dirname="blabla"

# assumes eos mounted on t3home
dirname="V20_emu_vs_V25"

rm -rf /eos/home-m/mratti/www/BHNL/recoAnalysis/signalRw/$dirname
cp -r plots/$dirname /eos/home-m/mratti/www/BHNL/recoAnalysis/signalRw/$dirname
cp HTACCESS /eos/home-m/mratti/www/BHNL/recoAnalysis/signalRw/$dirname/.htaccess

