CFLAGS = -pthread -m64 -Wno-deprecated -I src/
CFLAGS += $(shell root-config --cflags --libs)

all : TrackletAna Tracklet_RecoVtx Tracklet_GenHadron

TrackletAna : src/TrackletAna.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/TrackletAna.cxx src/Hit.h src/Tracklet.h -o TrackletAna

Tracklet_RecoVtx : src/Tracklet_RecoVtx.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/Tracklet_RecoVtx.cxx src/Hit.h src/Tracklet.h -o Tracklet_RecoVtx

Tracklet_GenHadron : src/Tracklet_GenHadron.cxx 
	g++ ${CFLAGS} src/Tracklet_GenHadron.cxx -o Tracklet_GenHadron