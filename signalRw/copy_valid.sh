dirname="1708_central_mainGrid"

# assumes eos mounted on t3home
dirname="1708_central_mainGrid"

rm -rf /eos/home-m/mratti/www/BHNL/recoAnalysis/signalValid/$dirname
cp -r plots/$dirname /eos/home-m/mratti/www/BHNL/recoAnalysis/signalValid/$dirname
cp HTACCESS /eos/home-m/mratti/www/BHNL/recoAnalysis/signalValid/$dirname/.htaccess

