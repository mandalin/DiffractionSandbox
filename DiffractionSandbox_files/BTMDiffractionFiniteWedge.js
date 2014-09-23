

function find_apex(rs, rr, srcz, recz){
var ClosestEdgePoint2Ppos=[0.0,0.0,recz];
var ClosestEdgePoint2Qpos=[0.0,0.0,srcz];
var zvec=subtract(ClosestEdgePoint2Ppos,ClosestEdgePoint2Qpos);
var Z=magnitude(zvec);
var zvecnorm=[zvec[0] / Z , zvec[1] / Z, zvec[2] / Z];
var za=Z*rs/(rs+rr); 
//from svensson2006, after nothing that zs=0, and zr=z for our setup
var apex_point=[za*zvecnorm[0]+ClosestEdgePoint2Qpos[0],za*zvecnorm[1]+ClosestEdgePoint2Qpos[1],za*zvecnorm[2]+ClosestEdgePoint2Qpos[2]];
return apex_point;
 };
 
function dot(left, right){
	var result
	result=left[0]*right[0] + left[1]*right[1] + left[2]*right[2];
	return result
};


function calc_least_time(r,ro,Z,co) { 
	//Z here is the distance between the source and receiver position on along the z axis. 
	least_time=Math.sqrt(((r+ro)*(r+ro) + Z*Z))/co; 
	return least_time;
 };
 
 
	function magnitude(pos1) { 
	mag=Math.sqrt(pos1[1]*pos1[1] + pos1[2]*pos1[2] + pos1[0]*pos1[0]); 
	return mag;
 };
 
 
 function subtract(pos1,pos2) { 
	var N=3
	var result = [];
	for (var i = 0; i < N; i++) 
	{
		result[i] = pos1[i]-pos2[i];
	}
	return result;
 };
 
 
function sinh(arg) {
  return (Math.exp(arg) - Math.exp(-arg)) / 2.0;
};

function acosh(x) {
  return Math.log(x + Math.sqrt(x * x - 1.0));
  //For all x>=1, we have arcosh(x)=ln(x+x2-1)
};

function which_case(apex, edge_end1, edge_end2){
	
            //   A  e1_________e2     case 1
			if ( (apex[2]<=edge_end1) && (apex[2]<=edge_end2) ) {	var whichcase=1;}
			
            //      e1____A____e2     case 2   
			if (((apex[2]<=edge_end1) && (apex[2]>=edge_end2))|| ((apex[2]>=edge_end1) && (apex[2]<=edge_end2))){var whichcase=2;}
				
            //      e1_________e2 A   case 3
			if ((apex[2]<=edge_end1) && (apex[2]<=edge_end2)) {var whichcase=3;}
			
			return whichcase;
};
 

 
 
 
 
//  //function BTM_IR( which_case,  theta_w, ro, thetao, r, theta,  Z, Ppos,  Qpos,  corner1,  corner2, least_time_point,  co,  rho,  S,  fs){		
//  function BTM_IR( which_case, corner1, corner2, Qpos, Ppos,co,least_time_point, theta_w,S,rho){		
//  	var fsdiff=4410.0;
//  	var Ndiff=Math.round(fsdiff);//this could be better...????!!!!

// 	//create empty IR
// 	var ImpulseResponse=[];
// 	ImpulseResponse.length=Ndiff;
// 	ImpulseResponse[0]=0;

//  	var infinite_wedge=false;
// 	var least_time;
//     var edge_distance_a, edge_time_a;
//     var edge_distance_b, edge_time_b;
//     var temp1, temp2;
// 	var edge_time_a, edge_time_b;
	
//     temp1=magnitude(subtract(corner1, Ppos));
//     temp2=magnitude(subtract(corner1, Qpos));
//     edge_distance_a = temp1 + temp2;
    
//  	temp1=magnitude(subtract(corner2, Ppos));
//     temp2=magnitude(subtract(corner2, Qpos));
//     edge_distance_b = temp1 + temp2;
    
//  	edge_time_a=edge_distance_a/co;
//     edge_time_b=edge_distance_b/co;

// 	var r=Math.sqrt(Ppos[0]*Ppos[0]+Ppos[1]*Ppos[0]);
// 	var ro=Math.sqrt(Qpos[0]*Qpos[0]+Qpos[1]*Qpos[0]);
// 	var Z=Math.sqrt((Ppos[2]-Qpos[2])*(Ppos[2]-Qpos[2]));
	
//  	least_time=calc_least_time(r,ro,Z,co);

//     var incident_sound_at_receiver_time=magnitude(subtract(Ppos,Qpos))/co;
//     var diffraction_delay=least_time-incident_sound_at_receiver_time;
    
    


// //     var t, tau, tauo, yarg, yargtau;
// //     var source_to_apex_distance=magnitude(subtract(Qpos,least_time_point));
// // 	var pi=Math.PI;
	

// // 	var thetao=Math.atan(Qpos[1]/Qpos[2]);
// // 	var theta=Math.atan(Ppos[1]/Ppos[2]);
	

	
// //     for(var n_diff=1;  n_diff<Ndiff ;n_diff++ ){ 
// // 		t=least_time+n_diff/fsdiff;  
// //  		tauo=least_time;
// //  		tau=t-tauo;
        
// //         yarg=((co*co*t*t) - ((r*r) + (ro*ro) + (Z*Z))) / (2.0*r*ro);
// // 		yargtau=((co*co*(2.0*(least_time*tau + tau*tau)))/(2.0*r*ro))+1;
        
// //          if(yarg  < 1){    
// //  			 ImpulseResponse[n_diff]=0;
// //          }
        
// //          else{
// //              y=acosh(yarg); 
            
            
// //              var Beta_num_1 = Math.sin( (pi / theta_w) * (pi + theta + thetao) );
// //              var Beta_num_2 = Math.sin( (pi / theta_w) * (pi + theta - thetao) );
// //              var Beta_num_3 = Math.sin( (pi / theta_w) * (pi - theta + thetao) );
// //              var Beta_num_4 = Math.sin( (pi / theta_w) * (pi - theta - thetao) );
            
// //  			 var exparg=(-pi * y) / theta_w;
// //  			 var cosarg=pi/theta_w;
			
// //              var Beta_denom_1=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi + theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
// //              var Beta_denom_2=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi + theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
// //              var Beta_denom_3=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi - theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);             Beta_denom_4=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi - theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w); 
            
// //              var Beta= (Beta_num_1/Beta_denom_1) + (Beta_num_2/Beta_denom_2) + (Beta_num_3/Beta_denom_3) + (Beta_num_4/Beta_denom_4);

            
            
// //              //for spherical waves
// //  			var without_approx=(r*ro*sinh(y));
// //  			var with_approx=Math.sqrt(2*tauo*tau*co*co*ro*r);
			
// //  			if(without_approx==0){
// //  	            p_of_t=((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/with_approx) * Math.exp((-pi*y) / theta_w);
// //  				}
// //  			else{
// //  	           	p_of_t=((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/without_approx) * Math.exp((-pi*y) / theta_w);
// //  		        }
				
// // //             //for plane waves 
// // //             //p_of_t=(source_to_apex_distance*4.0*pi) *((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/(r*ro*std::sinh(y))) * Math.exp((-pi*y) / theta_w);
            
// //              p_of_t=p_of_t/fsdiff;            

//              //   A  e1_________e2     case 1
//              //      e1____A____e2     case 2   
//              //      e1_________e2 A   case 3
            
// //             if(! infinite_wedge){
// //                 switch (which_case){
// //                     case 1:
// //                         if((t<edge_time_a)||(t>edge_time_b))  
// //                         {   ImpulseResponse[n_diff]=0;}
// //                         else 
// //                         {   ImpulseResponse[n_diff]=p_of_t/2;}
// //                         break;
                        
// //                     case 2:
// //                         if((t>edge_time_a) || (t>edge_time_b))
// //                         {   ImpulseResponse[n_diff]=p_of_t/2;}
// //                         else if((t>edge_time_a) && (t>edge_time_b))
// //                         {   ImpulseResponse[n_diff]=0;}
// //                         else
// //                         {   ImpulseResponse[n_diff]=p_of_t;}
// //                         break;
                        
// //                     case 3:
// //                         if((t<edge_time_b)||(t>edge_time_a)) 
// //                         {   ImpulseResponse[n_diff]=0;}
// //                         else
// //                         {   ImpulseResponse[n_diff]=p_of_t/2;}
// //                         break;
                        
// //                     default: cout << "UnspecifiedCase";
// //                         break;
// //                 }
// //             }
            
//          }
        
//      }
	
		
		
		
// 	return ImpulseResponse;
		
		
// };


