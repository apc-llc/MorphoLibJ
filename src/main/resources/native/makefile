ifeq (,$(JAVA_HOME))
  $(error "$$JAVA_HOME must be pointing to the JDK, e.g. $$ export JAVA_HOME=/opt/jdk1.8.0_102/)
endif

all: libWatershed.so

libWatershed.so: applyWithMask.cpp
	g++ -g -O3 -ffast-math -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux $< -fPIC -shared -o $@

clean:
	rm -rf libWatershed.so

