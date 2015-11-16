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

	void start(uint8_t p);

	void stop();

	void setStepSize(uint8_t s);

	uint8_t getStepSize();

	void update(double dir);




private:
	double pos_d;
	uint8_t pos;
	uint8_t step_size;

	uint8_t move_servo_pos[12];
	uint8_t lift_servo_pos[3];

	Servo move;
	Servo lift;

	bool started;


	


	/* data */
};





/* lengths between joints */
#define R1 0.05
#define R2 0.05
#define R3 0.05


/*
*
*								
*								/
*							   /\
*							  / psi
*							(C)---\---------(D)
*							/			R3
*						 R2/
*						  /
*			R1			 / \
*	(A)----------------(B)_ phi _ 
*	theta
*
*
*
*
*
*
*
*/






class RobotLeg3
{
public:
	RobotLeg3();

	void init(int s1_pin, int s2_pin);

	void start(uint8_t p);

	void stop();

	void setStepSize(uint8_t s);

	uint8_t getStepSize();

	void update(double dir);




private:
	/* help variables */
	matrix *Rx_theta, *Ry_phi, *Ry_psi;
	matrix *movement, *to_tip, *to_target;
	matrix *Z;
	matrix *R_Z, *R_R_Z, *R_R_R_Z;

	/* positions */
	matrix* A, B, C, tip, target;

	/* joint rotation axises */
	matrix* A_axis, B_axis, C_axis;

	/* Joint angles*/
	double theta, phi, psi;

	double gradient;




	double pos_d;
	uint8_t pos;
	uint8_t step_size;

	uint8_t move_servo_pos[12];
	uint8_t lift_servo_pos[3];

	Servo move;
	Servo lift;

	bool started;


	


	/* data */
};