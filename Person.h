class Person {
    
    // all members are public
public:
    
    int age;
    int gender;
    int stage;
    int zombieLife;
    int breedingCounter;

    int owner;
    int getOwner(){return owner;}
    void setOwner(int o){owner = o;}

    bool empty;
    bool isEmpty(){return empty;}
    void setEmpty(bool e){empty = e;}

    Person(){
        zombieLife = ZOMBIE_LIFE;
        breedingCounter = 0;
    }
    
    // convenience method
    bool canBreed(){
        return stage == HUMAN && empty == false && (age == ADULT || age == ELDERLY);
    }
    
    void decrementBreedingCounter(){
        if (breedingCounter > 0){
            breedingCounter -= 1;
        }
    }
    
    void startBreedingCounter(){
        breedingCounter = 270;
    }
    
    void setAge(int a){
        age = a;
    }
    
    int getAge(){
        return age;
    }
    
    void setGender(int g){
        gender = g;
    }
    
    int getGender(){
        return gender;
    }

    // EPIDEMIC ATTRIBUTES
    // 0 = Human, 1 = Exposed, 2 = Zombie, 3 = Dead
    
    // sets the stage
    void setStage(int s){
        stage = s;

        if (stage == ZOMBIE){
            setZombieLife(ZOMBIE_LIFE);
        }

    }
    
    // gets the stage
    int getStage(){
        return stage;
    }
    
    // advances the stage
    void advanceStage(){
        stage += 1;
    }
    
    // ZOMBIE ATTRIBUTES
    
    // sets the zombie life
    void setZombieLife(int days){
        zombieLife = days;
    }
    
    // gets the zombie life
    int getZombieLife(){
        return zombieLife;
    }
    
    // decreases the zombie life by one
    void decreaseZombieLife(){
        
        zombieLife -= 1;
        
        if (zombieLife <= 0){
            stage = DEAD;
            empty = true;
        }
        
    }

    // convenience method
    bool isZombie(){
        return stage == ZOMBIE;
    }

    bool isHuman(){
        return stage == HUMAN;
    }
    
    double getMoveAdjustment(double baseProb){
        
        // young move a little slower
        if (stage == 0){
            return baseProb - 0.025;
            
        // old move a lot slower
        } else if (stage == 2){
            return baseProb - 0.05;
            
        } else {
            return baseProb;
        }
        
    }

    double getInfectionProbability(){

        if (age == YOUNG){
            return PROB_INFECTION_YOUNG;
        } else if (age == ELDERLY){
            return PROB_INFECTION_ELDERLY;
        } else {
            return PROB_INFECTION;
        }

    }
};
