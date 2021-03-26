flag = -check all

objs = utility.o main.o

steqr.exe: $(objs)
	ifort -qopenmp -parallel -mkl -static-intel -ipo $(flag) $^ -o $@

%.o: source/%.f90
	ifort -qopenmp -parallel -mkl -static-intel -ipo $(flag) -c $<

.PHONY: clean
clean:
	rm -f *.mod *.o *.exe