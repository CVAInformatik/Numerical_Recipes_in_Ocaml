module TOEPLZ = struct

(*


	 This is based on Wikipedias description of Levinson Recursion and not the 'fortran-in-c syntax' version
	 from "Numerical Recipes in C, 2nd.edition" (Press et al)
	 
	 I first made an implemtation in C++ and when that worked, I ported it to Ocaml.
	 
	 Not too happy with my imperative programming skills in Ocaml...
 *)

exception SingularPrincipalMinor

let rec toeplz__  (r:float array) (x:float array) (y:float array) (n:int) (f:float array) (b:float array) inx = 
  if( inx = n) then ()
  else 
			begin
			  let enf = ref 0.0 in
			  let enb = ref 0.0 in
			  for t = 0 to (inx-1) do enf := !enf +. ( r.(n-1+inx-t) *. f.(t)) done; 
			  for t = 0 to (inx-1) do enb := !enb +. ( r.(n-2-t) *. b.(t)) done; 
			  let divisor = 1. -. (!enf *. !enb)
			  in if (( Float.classify_float divisor) <> Float.FP_normal) 
			  		then
			          raise SingularPrincipalMinor
			  		else
			  		 begin
			  		     for t = 0 to (inx-1) do f.(t) <- f.(t) /. divisor done; 
			  		     f.(inx) <- 0.0;
			  		     for t = inx downto 1 do b.(t) <- b.(t-1) /. divisor done; 
			  		     b.(0) <- 0.0;
			  		     for tt = 0 to inx do
			  		        let newf = f.(tt) -. (!enf *. b.(tt)) in
			  		        let newb = b.(tt) -. (!enb *. f.(tt)) in
			  		        f.(tt) <- newf;
			  		        b.(tt) <- newb;
			  		     done;
			  		     let enx = ref 0.0 in
			  		     begin
			  		     	for xx = 0 to (inx-1) do enx := (!enx) +. (r.(n-1+inx-xx) *. x.(xx)) done;
			  		     	for yy = 0 to (inx-1) do x.(yy)<- x.(yy) +. (b.(yy) *. ( y.(inx) -. !enx) ) done;
			  		      x.(inx) <- ((y.(inx) -. !enx) *. b.(inx));
									toeplz__ r x y n f b (inx+1);
			  		     end
			  		 end
				end			

(*
		Alternative version to toeplz if you know N > 2 and you have no problem with having to allocate
		the f and b arrays yourself
*)

 
let rec toeplz_  (r:float array) (x:float array) (y:float array) (n:int) (f:float array) (b:float array) = 
					          let determinant = (r.(n-1) *. r.(n-1)) -.( r.(n) *. r.(n-2))
					          in begin 
					          		x.(0)<- (( y.(0) *. r.(n-1)) -. ( y.(1) *. r.(n-2))) /. determinant ;
					           		x.(1)<- (( r.(n-1) *. y.(1)) -. ( y.(0) *. r.(n))) /. determinant ;
					           		b.(0)<- ( -1. *.  r.(n-2)) /. determinant ;
					           		b.(1)<-  r.(n-1)  /. determinant ;
					           		f.(0)<- ( r.(n-1))  /. determinant ;
					           		f.(1)<- ( -1. *. r.(n)) /. determinant ;					             	
					           		toeplz__  r x y n f b 2 ;
					            end

				 
(*
			 the parameters are in the same order and format as for "Num. Rec. in C", which in not necessarily
			 optimal for Ocaml.....
			 
			 r the coefficients of the toepliz matrix, size 2n-1,  
			        tminus1 is stored in r[n-2], t0 is stored in r[n-1], t1 in [n], t2 in [n+1].....
			 x the solution size n
			 y the right hand side size n
			 n the size of the arrays
			 
*)
	
let rec toeplz  (r:float array) (x:float array) (y:float array) (n:int) = 
					match n with
					| 0 ->  ()
					| 1 ->  begin if(( Float.classify_float r.(0)) = Float.FP_normal) 
									then 
											 x.(0) <-  (y.(0)) /. (r.(0))
					        else raise SingularPrincipalMinor 
					        end
					| 2 ->  let determinant = (r.(n-1) *. r.(n-1)) -. (r.(n) *. r.(n-2)) in
					        if (( Float.classify_float determinant) = Float.FP_normal) then
					             begin
					             	x.(0)<- (( y.(0) *. r.(n-1)) -. ( y.(1) *. r.(n-2))) /. determinant ;
					             	x.(1)<- (( r.(n-1) *. y.(1)) -. ( y.(0) *. r.(n  ))) /. determinant 
					             end
					        else begin
					           raise SingularPrincipalMinor
					        end			        					        
					| _ ->	if ( n> 2) then toeplz_ r x y n (Array.make n 0.0 ) (Array.make n 0.0 )
	
	
end 	

