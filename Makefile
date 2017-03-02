default:
	cd geometryEngines/TSurf && $(MAKE)
	cd meshTools/hypsurf && $(MAKE)
	cd utilities/CGNSinterface && $(MAKE)

clean:
	cd geometryEngines/TSurf && $(MAKE) clean
	cd meshTools/hypsurf && $(MAKE) clean
	cd utilities/CGNSinterface && $(MAKE) clean
