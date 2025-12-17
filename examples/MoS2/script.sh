#!/bin/bash

run_step() {
    local label="$1"
    local outfile="$2"
    shift 2

    echo 
    echo "[$label] started at $(date)"
    local t0
    t0=$(date +%s)

    "$@" > "$outfile" 2>&1
    local status=$?

    local t1
    t1=$(date +%s)
    local dt=$((t1 - t0))
    local h=$((dt / 3600))
    local m=$(((dt % 3600) / 60))
    local s=$((dt % 60))

    if [ $status -ne 0 ]; then
        echo "ERROR in [$label]: exit code $status"
        echo "See '$outfile' for details"
        exit $status
    fi

    printf "[%s] finished in %02d:%02d:%02d\n" "$label" "$h" "$m" "$s"
}

cd 01-scf
run_step "SCF" scf.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -in scf.in
cd ..

cd 02-wfn
cp -r ../01-scf/MoS2.save .
kgrid.x kgrid.in kgrid.out kgrid.log
# cat kgrid.out >> wfn.in
run_step "WFN" wfn.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -pd .true. -in wfn.in
run_step "WFN2BGW" pw2bgw.out \
    mpirun -np 1 -ppn 1 pw2bgw.x -pd .true. -in pw2bgw.in
rm -rf MoS2.save
run_step "WFN2HDF" wfn2hdf.out \
    wfn2hdf.x BIN WFN WFN.h5
rm -f WFN
cd ..

cd 03-wfnq
cp -r ../01-scf/MoS2.save .
setup_subsampling_nns.x HDF5 ../02-wfn/WFN.h5 3 > setup_subsampling_nns.out
# cat kpoints_all.dat >> wfn.in
run_step "WFNq" wfn.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -pd .true. -in wfn.in
run_step "WFNq2BGW" pw2bgw.out \
    mpirun -np 1 -ppn 1 pw2bgw.x -pd .true. -in pw2bgw.in
rm -rf MoS2.save
run_step "WFNq2HDF" wfn2hdf.out \
    wfn2hdf.x BIN WFN WFN.h5
rm -f WFN
cd ..

cd 04-wfn_co
cp -r ../01-scf/MoS2.save .
kgrid.x kgrid.in kgrid.out kgrid.log
# cat kgrid.out >> wfn.in
run_step "WFN_co" wfn.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -pd .true. -in wfn.in
run_step "WFN_co2BGW" pw2bgw.out \
    mpirun -np 1 -ppn 1 pw2bgw.x -pd .true. -in pw2bgw.in
rm -rf MoS2.save
run_step "WFN_co2HDF" wfn2hdf.out \
    wfn2hdf.x BIN WFN WFN.h5
rm -f WFN
cd ..

cd 05-wfnq_co
cp -r ../01-scf/MoS2.save .
kgrid.x kgrid.in kgrid.out kgrid.log
# cat kgrid.out >> wfn.in
run_step "WFNq_co" wfn.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -pd .true. -in wfn.in
run_step "WFNq_co2BGW" pw2bgw.out \
    mpirun -np 1 -ppn 1 pw2bgw.x -pd .true. -in pw2bgw.in
rm -rf MoS2.save
run_step "WFNq_co2HDF" wfn2hdf.out \
    wfn2hdf.x BIN WFN WFN.h5
rm -f WFN
cd ..

cd 06-wfn_fi
cp -r ../01-scf/MoS2.save .
kgrid.x kgrid.in kgrid.out kgrid.log
# echo 'K_POINTS crystal' >> wfn.in
# awk '/k-points in the original uniform grid/ {getline; print; exit}' kgrid.log >> wfn.in
# awk '/k-points in the original uniform grid/{f=1; next} f && !NF {f=0} f && ++c>1 {printf "%12.9f %12.9f %12.9f   1.0\n", $2, $3, $4}' kgrid.log >> wfn.in
run_step "WFN_fi" wfn.out \
    mpirun -np "$NP" -ppn "$PPN" pw.x -pd .true. -in wfn.in
run_step "wfn_fi2BGW" pw2bgw.out \
    mpirun -np 1 -ppn 1 pw2bgw.x -pd .true. -in pw2bgw.in
rm -rf MoS2.save
run_step "WFN_fi2HDF" wfn2hdf.out \
    wfn2hdf.x BIN WFN WFN.h5
rm -f WFN
cd ..

cd 07-epsilon
ln -sf ../02-wfn/WFN.h5
ln -sf ../03-wfnq/WFN.h5 WFNq.h5
# cat ../03-wfnq/epsilon_q0s.inp >> epsilon.inp
# awk 'NR>3 {printf " %12.9f %12.9f %12.9f   1.0 0\n", $1, $2, $3 }' ../02-wfn/kgrid.out | sed -i '/^end$/e cat' epsilon.inp
run_step "EPSILON" epsilon.out \
    mpirun -np "$NP" -ppn "$PPN" epsilon.cplx.x -in epsilon.inp
cd ..

cd 08-sigma
ln -sf ../02-wfn/WFN.h5 WFN_inner.h5
ln -sf ../02-wfn/RHO
ln -sf ../02-wfn/VXC
ln -sf ../03-wfnq/subweights.dat
ln -sf ../07-epsilon/eps0mat.h5
ln -sf ../07-epsilon/epsmat.h5
# awk 'NR>3 {printf " %12.9f %12.9f %12.9f   1.0\n", $1, $2, $3 }' ../02-wfn/kgrid.out >> sigma.inp
# echo 'end' >> sigma.inp
run_step "SIGMA" sigma.out \
    mpirun -np "$NP" -ppn "$PPN" sigma.cplx.x -in sigma.inp
cd ..

cd 09-epsilon_co
ln -sf ../04-wfn_co/WFN.h5
ln -sf ../05-wfnq_co/WFN.h5 WFNq.h5
# awk 'NR>3 {printf " %12.9f %12.9f %12.9f   1.0 0\n", $1, $2, $3 }' ../04-wfn_co/kgrid.out >> epsilon.inp
# echo 'end' >> epsilon.inp
run_step "EPSILON_co" epsilon.out \
    mpirun -np "$NP" -ppn "$PPN" epsilon.cplx.x -in epsilon.inp
cd ..

cd 10-inteqp
ln -sf ../02-wfn/WFN.h5 WFN_co.h5
ln -sf ../04-wfn_co/WFN.h5 WFN_fi.h5
ln -sf ../08-sigma/eqp1.dat eqp_co.dat
run_step "INTEQP" kernel.out \
    mpirun -np "$NP" -ppn "$PPN" inteqp.cplx.x -in inteqp.inp
cd ..

cd 11-kernel
ln -sf ../04-wfn_co/WFN.h5 WFN_co.h5
ln -sf ../09-epsilon_co/eps0mat.h5
ln -sf ../09-epsilon_co/epsmat.h5
run_step "KERNEL" kernel.out \
    mpirun -np "$NP" -ppn "$PPN" kernel.cplx.x -in kernel.inp
cd ..

cd 12-absorption
ln -sf ../04-wfn_co/WFN.h5 WFN_co.h5
ln -sf ../06-wfn_fi/WFN.h5 WFN_fi.h5
ln -sf ../10-inteqp/eqp.dat eqp_co.dat
ln -sf ../11-kernel/eps0mat.h5
ln -sf ../11-kernel/epsmat.h5
ln -sf ../11-kernel/bsemat.h5
run_step "ABSORPTION" absorption.out \
    mpirun -np "$NP" -ppn "$PPN" absorption.cplx.x -in absorption.inp
cd ..

