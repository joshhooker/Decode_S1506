root -l<<EOF
gSystem->Load("/home/reactions/treeIris/lib/libTEvent");
.L decode.cxx
decode($1)
EOF
