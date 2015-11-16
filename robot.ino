/* 
*/

#include "robotleg.h"

RobotLeg2 leg1; 

int potpin = 0; 
int val;   

void setup()
{
  leg1.init(8,9);
  Serial.begin(38400);
  delay(2000);

  leg1.start(6);

  delay(2);
}

void loop() 
{
	val = analogRead(potpin);
	val = map(val, 0, 1023, -10, 10);
	Serial.println(val);

	leg1.update((double)val/10.0);

	delay(20);

	

  
} 
