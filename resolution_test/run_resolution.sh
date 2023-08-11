echo "Original ident:"
read ident_org
echo "New ident: (e.g. resolutiontest)"
read ident_new

ln -s mat.$ident_org mat.$ident_new
ln -s ata.1.$ident_org ata.1.$ident_new
ln -s atamx.$ident_org atamx.$ident_new
ln -s columndensity.1.$ident_org columndensity.1.$ident_new
ln -s aux.$ident_org aux.$ident_new
