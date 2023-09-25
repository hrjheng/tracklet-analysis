CFLAGS = -pthread -m64 -Wno-deprecated -I src/
CFLAGS += $(shell root-config --cflags --libs)

all : Tracklet_BeamSpot Tracklet_BeamSpot_v2 TrackletAna Tracklet_RecoVtx plotRecoCluster plotTracklets Corrections ClusterRecoIneff

Tracklet_BeamSpot : src/Tracklet_BeamSpot.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/Tracklet_BeamSpot.cxx src/Hit.h src/Tracklet.h -o Tracklet_BeamSpot

Tracklet_BeamSpot_v2 : src/Tracklet_BeamSpot_v2.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/Tracklet_BeamSpot_v2.cxx src/Hit.h src/Tracklet.h -o Tracklet_BeamSpot_v2

TrackletAna : src/TrackletAna.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/TrackletAna.cxx src/Hit.h src/Tracklet.h -o TrackletAna

Tracklet_RecoVtx : src/Tracklet_RecoVtx.cxx src/Hit.h src/Tracklet.h
	g++ ${CFLAGS} src/Tracklet_RecoVtx.cxx src/Hit.h src/Tracklet.h -o Tracklet_RecoVtx

plotRecoCluster : src/plotRecoCluster.cxx
	g++ $(CFLAGS) src/plotRecoCluster.cxx src/Hit.h src/Tracklet.h src/GenHadron.h -o plotRecoCluster

plotTracklets : src/plotTracklets.cxx
	g++ $(CFLAGS) src/plotTracklets.cxx src/Hit.h src/Tracklet.h src/GenHadron.h -o plotTracklets

Corrections : src/Corrections.cxx
	g++ $(CFLAGS) src/Corrections.cxx -o Corrections

ClusterRecoIneff : src/ClusterRecoIneff.cxx
	g++ $(CFLAGS) src/ClusterRecoIneff.cxx -o ClusterRecoIneff