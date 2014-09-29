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
			if ( (apex[2]<=edge_end1[2]) && (apex[2]<=edge_end2[2]) ) {	var whichcase=1;}
			
            //      e1____A____e2     case 2   
			if (((apex[2]<=edge_end1[2]) && (apex[2]>=edge_end2[2]))|| ((apex[2]>=edge_end1[2]) && (apex[2]<=edge_end2[2]))){var whichcase=2;}
				
            //      e1_________e2 A   case 3
			if ((apex[2]<=edge_end1[2]) && (apex[2]<=edge_end2[2])) {var whichcase=3;}
			
			return whichcase;
};
 
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
 
 
 function BTM_IR( which_case, Ppos, Qpos, corner_1, corner_2, co,least_time, fsdiff, apex, rr, rs, Z, theta_w){		
	var infinite_wedge=false;
    var least_time;
    var edge_distance_a, edge_time_a;
    var edge_distance_b, edge_time_b;
    var temp1, temp2;
	var edge_time_a, edge_time_b;
	

    temp1=magnitude(subtract(corner_1, Ppos));
    temp2=magnitude(subtract(corner_1, Qpos));
    edge_distance_a = temp1 + temp2;
   	
   	temp1=magnitude(subtract(corner_2, Ppos));
    temp2=magnitude(subtract(corner_2, Qpos));
	edge_distance_b = temp1 + temp2;
    
 	edge_time_a=edge_distance_a/co;
	edge_time_b=edge_distance_b/co;
	
	var incident_sound_at_receiver_time=magnitude(subtract(Ppos,Qpos))/co;
    var diffraction_delay=least_time-incident_sound_at_receiver_time;
	var diffraction_delay_samples=diffraction_delay/fsdiff;
   

    var yarg, y, Beta_num_1, Beta_num_2, Beta_num_3, Beta_num_4 , Beta_denom_1, Beta_denom_2, Beta_denom_3, Beta_denom_4, Beta, p_of_t;
    
    
  	//Create Empty IR
    var Ndiff=Math.round(fsdiff/100);//this could be better...????!!!!
   	//Ndiff=100; //For Fig8 Medwin 1982
	var ImpulseResponse=[];
	ImpulseResponse.length=Ndiff;

 
  	var longest_time;
    if(edge_time_a>edge_time_b){   
 		longest_time=edge_time_a;
 	}else{ longest_time=edge_time_b;
 	}

    
    var t, tau, tauo, yargtau, theta, thetao;
    var source_to_apex_distance=magnitude(subtract(Qpos,apex));
 	var pi=Math.PI;


	var AlongFaceUnitVector=[-1,0,0];
	
 	var edge_2_Ppos=[Ppos[0], Ppos[1], 0];
	var edge_2_Ppos_unit=[edge_2_Ppos[0]/magnitude(edge_2_Ppos), edge_2_Ppos[1]/magnitude(edge_2_Ppos),0];
	var theta_temp=dot(AlongFaceUnitVector,edge_2_Ppos_unit);
	var theta_acute=Math.acos(theta_temp);
	if(Ppos[1]<0){ theta=2*pi-theta_acute;}
	else{theta=theta_acute;}
	
 	var edge_2_Qpos=[Qpos[0], Qpos[1], 0];
 	var edge_2_Qpos_unit=[edge_2_Qpos[0]/magnitude(edge_2_Qpos), edge_2_Qpos[1]/magnitude(edge_2_Qpos),0];
 	var thetao_temp=dot(AlongFaceUnitVector,edge_2_Qpos_unit);
 	var thetao_acute=Math.acos(thetao_temp);
	if(Qpos[1]<0){ thetao=2*pi-thetao_acute;}
	else{thetao=thetao_acute;}
	
			
	

	
 	for(var n_diff=0;  n_diff<Ndiff ;n_diff++ ){ 
  		t=least_time+n_diff/fsdiff+(.5/fsdiff);  
		tau=t-least_time;
		tauo=least_time;
		
        
		yarg=((co*co*t*t) - ((rr*rr) + (rs*rs) + (Z*Z))) / (2.0*rr*rs);
		yargtau=((co*co*(2.0*(least_time*tau + tau*tau)))/(2.0*rr*rs))+1;
        
		if(yarg  < 1){    
			ImpulseResponse[n_diff]=0;
		}
        
		else{
         	y=acosh(yarg); 

			var Beta_num_1 = Math.sin( (pi / theta_w) * (pi + theta + thetao) );
			var Beta_num_2 = Math.sin( (pi / theta_w) * (pi + theta - thetao) );
			var Beta_num_3 = Math.sin( (pi / theta_w) * (pi - theta + thetao) );
			var Beta_num_4 = Math.sin( (pi / theta_w) * (pi - theta - thetao) );            
		}
		
 			var exparg=(-pi * y) / theta_w;
 			var cosarg=pi/theta_w;
			
             var Beta_denom_1=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi + theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
             var Beta_denom_2=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi + theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
             var Beta_denom_3=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi - theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
             var Beta_denom_4=1.0-(2.0* Math.exp(exparg) * Math.cos( cosarg * (pi - theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w); 
            
             Beta= (Beta_num_1/Beta_denom_1) + (Beta_num_2/Beta_denom_2) + (Beta_num_3/Beta_denom_3) + (Beta_num_4/Beta_denom_4);

            
            
//             //for spherical waves
			 var without_approximation=(rr*rs*sinh(y));
			 var with_approximation=Math.sqrt(2*tauo*(tau)*co*co*rs*rr);
			 S=1;
			 rho=1.2;
			 
			 //if(without_approximation==0){  
			 //		p_of_t=((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/with_approximation) * Math.exp((-pi*y) / theta_w);
			 //}else{
				 
			if (n_diff==0){
				p_of_t=0;
			}else{
	             p_of_t=((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/without_approximation) * Math.exp((-pi*y) / theta_w);
			}
			 //}
			 
//             //for plane waves 
//             //p_of_t=(source_to_apex_distance*4.0*pi) *((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/(r*ro*std::sinh(y))) * Math.exp((-pi*y) / theta_w);
            
             p_of_t=p_of_t/fsdiff;            

            //   A  e1_________e2     case 1
            //      e1____A____e2     case 2   
            //      e1_________e2 A   case 3
            
            if(! infinite_wedge){
                switch (which_case){
                    case 1:
                        if((t<edge_time_a)||(t>edge_time_b))  
                        {   ImpulseResponse[n_diff]=0;}
                        else 
                        {   ImpulseResponse[n_diff]=p_of_t/2;}
                        break;
                        
                    case 2:
//                          if((t>edge_time_a) || (t>edge_time_b))
//                          {   ImpulseResponse[n_diff]=p_of_t/2;}
//                          else if((t>edge_time_a) && (t>edge_time_b))
//                          {   ImpulseResponse[n_diff]=0;}
//                          else
                        {   ImpulseResponse[n_diff]=p_of_t;}
                        break;
                        
                    case 3:
                        if((t<edge_time_b)||(t>edge_time_a)) 
                        {   ImpulseResponse[n_diff]=0;}
                        else
                        {   ImpulseResponse[n_diff]=p_of_t/2;}
                        break;
                        
                    default: cout << "UnspecifiedCase";
                        break;
                }
            }
        
   }
	
		
	return ImpulseResponse;
		
		
};





    

///////////////

	
//
//     
//							   //   position_vector PposDiff, QposDiff;
//								//  PposDiff.assign(12.2,12.21,.01);    
//							   //   QposDiff.assign(12.2,12.21,.01);  
//    
//    
//    var infinite_wedge=false;
//    var least_time;
//    var edge_distance_a, edge_time_a;
//    var edge_distance_b, edge_time_b;
//    var temp1, temp2;
//    temp1=magnitude(sub(corner1, Ppos));
//    temp2=magnitude(sub(corner1, Qpos));
//    edge_distance_a = temp1 + temp2;
//    temp1=magnitude(sub(corner2, Ppos));
//    temp2=magnitude(sub(corner2, Qpos));
//    edge_distance_b = temp1 + temp2;
//    edge_time_a=edge_distance_a/co;
//    edge_time_b=edge_distance_b/co;
//    
//								////////////////////////////////////////////////////////    
//								//      This is For Validation      Fig 2 Medwin 1982 //
//								////////////////////////////////////////////////////////   
//								//////////////////////////////////////////////////////////////////////
//								//    S=1.0;
//								//    rho=1.2;                                                        //
//								//    Which_IR__=1;                                                   //
//								//    theta_w=(pi/180.00)*270.0;                                      //
//								//    // double ro, thetao, r, theta, Z;                              //
//								//    //Source Coordinant in cylindrical (ro, thetao, 0) from Qpos    //
//								//    ro=100.0;                                                       //
//								//    thetao=pi/4.0;                                                  //  
//								//    //Receiver Coordinant in cylindrical (r, theta, Z) from Ppos    //
//								//    double epsilon=(pi/180.0)*.1;                                   //  
//								//    theta=((pi/180.0)*225.0)+epsilon;                               //
//								//    r=100.0;                                                        //
//								//    Z=0;                                                            //
//								//    co=340;                                                         //
//								//    fsdiff=100000;                                                   //
//								//    infinite_wedge=true;                                            //
//								//////////////////////////////////////////////////////////////////////
//									
//									
//								////////////////////////////////////////////////////////    
//								//      This is For Validation      Fig 8 Medwin 1982 //
//								////////////////////////////////////////////////////////   
//								//////////////////////////////////////////////////////////////////////
//								//    S=1.0;                                                          //
//								//    rho=1.2;                                                        //
//								//    Which_IR__=1;                                                   //
//								//    theta_w=(pi/180.00)*360.0;                                      //
//								//    // double ro, thetao, r, theta, Z;                              //
//								//    //Source Coordinant in cylindrical (ro, thetao, 0) from Qpos    //
//								//    ro=.25;                                                         //
//								//    thetao=(pi/180.0)*15.0;                                         //  
//								//    //Receiver Coordinant in cylindrical (r, theta, Z) from Ppos    //                                    
//								//    theta=((pi/180.0)*225.0);                                       //
//								//    r=.25;                                                          //
//								//    Z=0;                                                            //
//								//    co=340;                                                         //
//								//    fsdiff=80000;                                                   //
//								//    infinite_wedge=true;                                            //
//								//////////////////////////////////////////////////////////////////////    
//    
//    least_time=(std::sqrt((r+ro)*(r+ro) + Z*Z))/co;
//
//    double incident_sound_at_receiver_time=magnitude(sub(Ppos,Qpos))/co;
//    double diffraction_delay=least_time-incident_sound_at_receiver_time;
//    
//
//    std::cout<<"For Simulation Number "<<Which_SIM__<<" Impulse Response Number "<< Which_IR__ <<std::endl;
//    std::cout<<"Also know as Simulation: "<<simulation_name<<std::endl;
//    std::cout<<"The diffraction delay is: "<<diffraction_delay;
//    
//    bool IR_written;
//    double yarg, y, Beta_num_1, Beta_num_2, Beta_num_3, Beta_num_4 , Beta_denom_1, Beta_denom_2, Beta_denom_3, Beta_denom_4, Beta, p_of_t;
//    
//    
//    //Create Empty IR
//    // int Ndiff=end_of_edge_time*fsdiff+10;  //good for finite edges
//    double longest_time;
//    if(edge_time_a>edge_time_b)
//    {   longest_time=edge_time_a;}
//    else
//    { longest_time=edge_time_b;}
//
//    int Ndiff=round(fsdiff);//this could be better...????!!!!
//    //Ndiff=100; //For Fig8 Medwin 1982
//    
//    std::vector<double> ImpulseResponse(Ndiff);
//    
//    double t, tau, yargtau;
//    double source_to_apex_distance=magnitude(sub(Qpos,least_time_point));
//    std::cout<<endl<<endl;
//    for(int n_diff=1;  n_diff<Ndiff ;n_diff++ )
//    { 
//        t=least_time+n_diff/fsdiff;  
//        
//        yarg=((co*co*t*t) - ((r*r) + (ro*ro) + (Z*Z))) / (2.0*r*ro);
//        //yargtau=((co*co*(2.0*(least_time*tau + tau*tau)))/(2.0*r*ro))+1;
//        
//        if(yarg  < 1)
//        {    std::cout<<"Dude, there is an error, because the argument for y is negative and we are not prepared for imaginary number results"<<endl;
//            ImpulseResponse[n_diff]=0;
//            // return IR_written;
//        }
//        
//        else
//        {
//            y=acosh(yarg);
//            
//            Beta_num_1 = sin( (pi / theta_w) * (pi + theta + thetao) );
//            Beta_num_2 = sin( (pi / theta_w) * (pi + theta - thetao) );
//            Beta_num_3 = sin( (pi / theta_w) * (pi - theta + thetao) );
//            Beta_num_4 = sin( (pi / theta_w) * (pi - theta - thetao) );
//            
//            Beta_denom_1=1.0-(2.0* Math.exp((-pi * y) / theta_w) * Math.cos( (pi/theta_w) * (pi + theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
//            Beta_denom_2=1.0-(2.0* Math.exp((-pi * y) / theta_w) * Math.cos( (pi/theta_w) * (pi + theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
//            Beta_denom_3=1.0-(2.0* Math.exp((-pi * y) / theta_w) * Math.cos( (pi/theta_w) * (pi - theta + thetao)  )) + Math.exp((-2.0*pi*y) / theta_w);
//            Beta_denom_4=1.0-(2.0* Math.exp((-pi * y) / theta_w) * Math.cos( (pi/theta_w) * (pi - theta - thetao)  )) + Math.exp((-2.0*pi*y) / theta_w); 
//            
//            Beta= (Beta_num_1/Beta_denom_1) + (Beta_num_2/Beta_denom_2) + (Beta_num_3/Beta_denom_3) + (Beta_num_4/Beta_denom_4);
//
//            
//            
//            //for spherical waves
//            //p_of_t=((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/(r*ro*std::sinh(y))) * Math.exp((-pi*y) / theta_w);
//
//            //for plane waves 
//            p_of_t=(source_to_apex_distance*4.0*pi) *((-S*rho*co)/(4.0*pi*theta_w)) * Beta * (1.0/(r*ro*std::sinh(y))) * Math.exp((-pi*y) / theta_w);
//            
//            p_of_t=p_of_t/fsdiff;            
//
//            //   A  e1_________e2     case 1
//            //      e1____A____e2     case 2   
//            //      e1_________e2 A   case 3
//            
//            if(! infinite_wedge)
//            {
//                switch (which_case)
//                {
//                    case 1:
//                        if((t<edge_time_a)||(t>edge_time_b))  
//                        {   ImpulseResponse[n_diff]=0;}
//                        else 
//                        {   ImpulseResponse[n_diff]=p_of_t/2;}
//                        break;
//                        
//                    case 2:
//                        if((t>edge_time_a) || (t>edge_time_b))
//                        {   ImpulseResponse[n_diff]=p_of_t/2;}
//                        else if((t>edge_time_a) && (t>edge_time_b))
//                        {   ImpulseResponse[n_diff]=0;}
//                        else
//                        {   ImpulseResponse[n_diff]=p_of_t;}
//                        break;
//                        
//                    case 3:
//                        if((t<edge_time_b)||(t>edge_time_a)) 
//                        {   ImpulseResponse[n_diff]=0;}
//                        else
//                        {   ImpulseResponse[n_diff]=p_of_t/2;}
//                        break;
//                        
//                    default: cout << "UnspecifiedCase";
//                        break;
//                }
//            }
//           std::cout<<p_of_t<<" ";
//            
//        }
//        
//    }
//    
//
//    //////////////Plane Wave Scaling/////////////////////////
//    double incident_pressure_at_apex,Prop_dist;
//    position_vector src_to_apex;
//    
//    src_to_apex=sub(least_time_point,Qpos);
//    Prop_dist=magnitude(src_to_apex);
//    
//    incident_pressure_at_apex=(rho*S/(4*pi*Prop_dist));
//    /////////////////////////////////////////////////////////
//    
//    double A,FirstSamp;
//    A= S*rho*Beta/(4*pi*theta_w*std::sqrt(2*least_time*ro*r));
//    FirstSamp=A/incident_pressure_at_apex;
//    
//
//    FILE * writeIR;
//    writeIR = fopen (simulation_name,"w");
//    IR_written=writeIR; //impulse response written successfully ?
//    
//    if(writeIR)
//    {
//        fprintf(writeIR,"%f \n",diffraction_delay  );
//        fprintf(writeIR,"%f \n\n",incident_pressure_at_apex  );
//        
//        //  fputs ("y=[",write);
//        for(unsigned sample=0; sample<Ndiff; sample++)
//        {
//            fprintf(writeIR,"%f \n",(ImpulseResponse[sample]));
//        }
//        // fprintf(write,"];");
//        
//        fclose(writeIR);
//        
//    }
//    else
//    {
//        std::cout<<"Couldnt Open FIle"<<std::endl;
//    }
//    
//    //----------------------------------------------------------------------------------------------------------
//    
//
//            
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    return IR_written;

    


//