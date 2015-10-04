

#include "leg2.h"






RobotLeg2::RobotLeg2(){
	


}

void RobotLeg2::init(int s1_pin, int s2_pin){


}

void RobotLeg2::setSetpSize(uint8_t s){
	step_size = s;

}

void RobotLeg2::update(uint8_t dir){
	pos = (pos + dir)%MAX_COUNT;

	if(pos < LIFT_PUT){ /* move phase */

		move.write(move_servo_pos[pos]);

	}else if(pos < RESET){ /* lift/put phase */
		switch(pos){
			case 12
				lift.write(lift_servo_pos[0]);
			break;
			case 13
				lift.write(lift_servo_pos[1]);
			break;
		}

	}else if(pos < PUT_LIFT){ /* reset phase */

		switch(pos){
			case 14
				move.write(move_servo_pos[11]);
			break;
			case 15
				move.write(move_servo_pos[7]);
			break;
			case 16
				move.write(move_servo_pos[3]);
			break;
			case 17
				move.write(move_servo_pos[0]);
			break;
		}


	}else{ /* put/lift phase */
		switch(pos){
			case 18
				lift.write(lift_servo_pos[1]);
			break;
			case 19
				lift.write(lift_servo_pos[0]);
			break;
		}		
	}


}