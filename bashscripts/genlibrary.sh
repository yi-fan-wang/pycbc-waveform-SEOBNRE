gcc -fPIC -shared src/main-interface.c src/SEOBNREforPythonMain.c -L/opt/local/lib -lgsl -lgslcblas -lm -o libSEOBNRE.so
