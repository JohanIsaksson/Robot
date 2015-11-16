

#include "robotleg.h"






RobotLeg2::RobotLeg2(){
	


}

void RobotLeg2::init(int s1_pin, int s2_pin){

	move.attach(s1_pin);
	lift.attach(s2_pin);
	stop();

	lift_servo_pos[0] = 70;
	lift_servo_pos[1] = 85;

	setStepSize(STANDARD_STEP);



}

/* requires som delay to let servos aquire correct position */
void RobotLeg2::start(uint8_t p){	
	started = true;
	pos = 6;
	pos_d = 6.0;
	update(0.0);
	
}


void RobotLeg2::stop(){
	started = false;
	move.write(90);
	lift.write(179);
}

void RobotLeg2::setStepSize(uint8_t s){
	step_size = s;
	for(int i = 0; i < 12; i++){
		move_servo_pos[i] = (90-s/2) + (uint8_t)((float)i * ((float)s/11.0));
	}

}

uint8_t RobotLeg2::getStepSize(){
	return step_size;
}

void RobotLeg2::update(double dir){

	if(pos_d + dir < 0.0){
		pos_d = 19.0;
	}else if(pos_d + dir > 19.0){
		pos_d = 0.0;
	}else{
		pos_d = pos_d + dir;
	}

	pos = (int8_t)pos_d;

	if(pos < LIFT_PUT){ /* move phase */

		move.write(move_servo_pos[pos]);
		lift.write(lift_servo_pos[0]);

	}else if(pos < RESET){ /* lift/put phase */
		switch(pos){
			case 12:
				lift.write(lift_servo_pos[0]);
			break;
			case 13:
				lift.write(lift_servo_pos[1]);
			break;
		}
		move.write(move_servo_pos[11]);

	}else if(pos < PUT_LIFT){ /* reset phase */

		switch(pos){
			case 14:
				move.write(move_servo_pos[11]);
			break;
			case 15:
				move.write(move_servo_pos[7]);
			break;
			case 16:
				move.write(move_servo_pos[3]);
			break;
			case 17:
				move.write(move_servo_pos[0]);
			break;
		}
		lift.write(lift_servo_pos[1]);


	}else{ /* put/lift phase */
		switch(pos){
			case 18:
				lift.write(lift_servo_pos[1]);
			break;
			case 19:
				lift.write(lift_servo_pos[0]);
			break;
		}		
		move.write(move_servo_pos[0]);
	}


}

void RobotLeg3::calculate_rotation_matrices(){
	value arr[9] = {
		cos(theta), 0, -sin(theta),
		0,			1,	0,
		sin(theta), 0, cos(theta)
	};
	insert_array(arr, Rx_theta);

	arr[9] = {
		1, 	0,			0,
		0, cos(phi), -sin(phi),
		0, sin(phi), cos(phi)
	};
	insert_array(arr, Ry_phi);

	arr[9] = {
		1, 	0,			0,
		0, cos(psi), -sin(psi),
		0, sin(psi), cos(psi)
	};
	insert_array(arr, Ry_psi);
}

void RobotLeg3::calculate_joint_positions(){

	/* precalculate rotations */
	calculate_rotation_matrices();

	multiply_matrices(Rx_theta, z, R_Z);
	multiply_matrices(Rx_phi, R_Z, R_R_Z);
	multiply_matrices(Rx_psi, R_R_Z, R_R_R_Z);

	/* B */
	multiply_matrices(Rx_theta, Z, R_Z);
	multiply_matrix_with_scalar(R1,R_Z);
	add_matrices(A, R_Z, B);

	/* C */
	multiply_matrices(Ry_phi, Z, R_Z);
	multiply_matrices(Rx_theta, R_Z, R_R_Z);
	multiply_matrix_with_scalar(R2, R_R_Z);
	add_matrices(B, R_R_Z, C);

	/* Tip */
	multiply_matrices(Ry_psi, Z, R_R_Z);
	multiply_matrices(Ry_phi, R_Z, R_R_Z);
	multiply_matrices(Rx_theta, R_R_Z, R_R_R_Z);
	multiply_matrix_with_scalar(R3, R_R_R_Z);
	add_matrices(C, R_R_R_Z, tip);
	
}

void RobotLeg3::gradient_follow(){

	subtract_matrices(target, tip, to_target);

	/* A */
		subtract_matrices(tip, A, to_tip);		

		cross_product(to_tip, A_axis, movement);
		gradient = dot_product(movement,to_target);


	/* B */

	/* C */

}



for each joint
        if 3D:  axis = axis of rotation for this joint 
        if 2D: 	axis = (0, 0, 1)

        ToTip = tip - joint_centre

        ToTarget = target - tip

        movement_vector = crossproduct(ToTip, axis)
        gradient = dotproduct(movement_vector, ToTarget)    
    end loop