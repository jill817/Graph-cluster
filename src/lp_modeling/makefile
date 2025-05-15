CC := g++ -std=c++11
sol := sls_solver/

nia_ls: nia_ls.o sls_solver/solver.o sls_solver/model.o
	$(CC) -O3 -o $@ $^

%.o: %.cpp
	$(CC) -c -O3 -o $@ $<

solve: nia_ls
	@echo ""
	@echo "###################################"
	@echo "# TEST: solve               #"
	@echo "###################################" 
	./nia_ls solve
	@echo ""
	@echo "###################################"
	@echo "# solve done!                   #"
	@echo "###################################"

clean:
	find . -name "*.o" | xargs rm -rf
	rm -f nia_ls

resolve: 
	make clean
	make solve

demo_solve: nia_ls
	@echo ""
	@echo "###################################"
	@echo "# TEST: demo solve               #"
	@echo "###################################" 
	./nia_ls dsolve
	@echo ""
	@echo "###################################"
	@echo "# demo solve done!                   #"
	@echo "###################################"

demo_lp: nia_ls
	@echo ""
	@echo "###################################"
	@echo "# TEST: demo write lp               #"
	@echo "###################################" 
	./nia_ls dlp
	@echo ""
	@echo "###################################"
	@echo "# demo write lp done!                   #"
	@echo "###################################"

clz:
	rm -rf experiment
	make clean

experiment: nia_ls
	mkdir experiment
	@echo ""
	@echo "###################################"
	@echo "# TEST: demo write lp               #"
	@echo "###################################" 
	sh run.sh
	@echo ""
	@echo "###################################"
	@echo "# demo write lp done!                   #"