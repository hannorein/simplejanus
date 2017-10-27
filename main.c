#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Scale for conversion between FP and INT
const double scale_vel  = 1e-16; 
const double scale_pos  = 1e-16;

// INT datatype used
#define REB_PARTICLE_INT_TYPE int64_t

struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};

double t = 0;
const double dt = 0.01;

const unsigned int N = 9;
struct reb_particle_int* p_int = NULL;

struct reb_particle {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    double m;
};

// Initial conditions (in FP)
struct reb_particle p[] = { 
    	{
		.x = 0.0021709922250528, .y = 0.0057845061154043, .z = -0.0001290326677066,
		.vx = -0.0003084904334499, .vy = 0.0003164862379414, .vz = 0.0000072860648107,
		.m = 1.0000000000000000
	},
	{
		.x = -0.1529074277548495, .y = -0.4329649810759809, .z = -0.0217536815870956,
		.vx = 1.2130892048755062, .vy = -0.4636664872138580, .vz = -0.1492230266727991,
		.m = 0.0000001660114153
	},
	{
		.x = -0.7051385792282048, .y = 0.1305062392874893, .z = 0.0423980407616931,
		.vx = -0.2107118903711193, .vy = -1.1628741220859935, .vz = -0.0038067721592922,
		.m = 0.0000024478382878
	},
	{
		.x = 0.8303864923760965, .y = 0.5551748865431479, .z = -0.0001556226179998,
		.vx = -0.5694403294004744, .vy = 0.8300359440285254, .vz = -0.0000250486216637,
		.m = 0.0000030404326480
	},
	{
		.x = -1.6007632981663540, .y = 0.4507843866326728, .z = 0.0485350310380760,
		.vx = -0.1874661855400607, .vy = -0.7140231189065021, .vz = -0.0103688562255236,
		.m = 0.0000003227156038
	},
	{
		.x = -4.5444724195553627, .y = -2.9811209359531872, .z = 0.1140115745580475,
		.vx = 0.2354668506120313, .vy = -0.3459544002171689, .vz = -0.0038305410200901,
		.m = 0.0009547919152112
	},
	{
		.x = -0.2998316596246585, .y = -10.0512228718170959, .z = 0.1866942196718307,
		.vx = 0.3063599906570191, .vy = -0.0107135147677418, .vz = -0.0120072161180579,
		.m = 0.0002858856727222
	},
	{
		.x = 17.8418531053445939, .y = 8.8433796310403689, .z = -0.1982994964737093,
		.vx = -0.1032131635550300, .vy = 0.1941992816066720, .vz = 0.0020584917278455,
		.m = 0.0000436624373583
	},
	{
		.x = 28.6228992820092181, .y = -8.7910334836014847, .z = -0.4786090163574258,
		.vx = 0.0523633993793736, .vy = 0.1755278382196959, .vz = -0.0048214129381180,
		.m = 0.0000515138377263
	},
};



static void to_int(){
    for(unsigned int i=0; i<N; i++){ 
        p_int[i].x = p[i].x/scale_pos; 
        p_int[i].y = p[i].y/scale_pos; 
        p_int[i].z = p[i].z/scale_pos; 
        p_int[i].vx = p[i].vx/scale_vel; 
        p_int[i].vy = p[i].vy/scale_vel; 
        p_int[i].vz = p[i].vz/scale_vel; 
    }
}
static void to_double(){
    for(unsigned int i=0; i<N; i++){ 
        p[i].x = ((double)p_int[i].x)*scale_pos; 
        p[i].y = ((double)p_int[i].y)*scale_pos; 
        p[i].z = ((double)p_int[i].z)*scale_pos; 
        p[i].vx = ((double)p_int[i].vx)*scale_vel; 
        p[i].vy = ((double)p_int[i].vy)*scale_vel; 
        p[i].vz = ((double)p_int[i].vz)*scale_vel; 
    }
}

static void drift(){
    for(unsigned int i=0; i<N; i++){
        p_int[i].x += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[i].vx*scale_vel/scale_pos);
        p_int[i].y += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[i].vy*scale_vel/scale_pos);
        p_int[i].z += (REB_PARTICLE_INT_TYPE)(dt/2.*(double)p_int[i].vz*scale_vel/scale_pos);
    }
}

static void kick(){
    for(unsigned int i=0; i<N; i++){
        p_int[i].vx += (REB_PARTICLE_INT_TYPE)(dt*p[i].ax/scale_vel);
        p_int[i].vy += (REB_PARTICLE_INT_TYPE)(dt*p[i].ay/scale_vel);
        p_int[i].vz += (REB_PARTICLE_INT_TYPE)(dt*p[i].az/scale_vel);
    }
}

static void gravity(){
    for(unsigned int i=0; i<N; i++){
        p[i].ax = 0.;
        p[i].ay = 0.;
        p[i].az = 0.;
        for(unsigned int j=0; j<N; j++){
            if (i!=j){
                const double dx = p[i].x - p[j].x;
                const double dy = p[i].y - p[j].y;
                const double dz = p[i].z - p[j].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz);
                const double prefact = -1/(_r*_r*_r)*p[j].m;
                
                p[i].ax    += prefact*dx;
                p[i].ay    += prefact*dy;
                p[i].az    += prefact*dz;
            }
        }
    }
}

void janus_step(){
    // One leapfrog step
    drift();
    
    to_double(); 
    gravity();
    kick();

    drift();

    t += dt;
    
    // Only needed for floating point outputs 
    to_double(); 
}

int main(){
    // Setup initial integer coordinates
    if (p_int==NULL){
        p_int = malloc(sizeof(struct reb_particle_int)*N);
        to_int(); 
    }
    
    while (t<2.*M_PI*1e1){ // 1 year
        janus_step();
        for(int i=0;i<N;i++){
            printf("%e %e\n",p[i].x,p[i].y);
        }
    }

    return 1;
}
