#include <Servo.h>


#define START 0
#define LIFT_PUT 12
#define RESET 14
#define PUT_LIFT 18

#define MAX_COUNT 20
#define STANDARD_STEP 40

class RobotLeg2
{
public:
	RobotLeg2();

	void init(int s1_pin, int s2_pin);

	void update();




private:
	uint8_t pos;
	uint8_t step_size;

	uint8_t move_servo_pos[12];
	uint8_t lift_servo_pos[2];

	Servo move;
	Servo lift;


	


	/* data */
};