/*
 * uart.xc
 *
 *  Created on: 2019.10.10
 *  Author: Xu Zongxiang
 */
#include "uart.h"
#include <xs1.h>
#include <platform.h>
#include <stdio.h>
#include <string.h>
#include "gpio_access.h"
#include <safe/string.h>
#include "flash_drive.h"
#include "protocol_frame.h"
#include "flash_drive.h"

#define BAUDRATE 115200
#define BIT_TIME XS1_TIMER_HZ / BAUDRATE
#define DATA_BITS 8

extern void device_reboot1(void);//rebot xmos

unsigned char wifi_id[13]="TA1911000101A";
unsigned char readdata[256]={0};
struct code_version{
    short BlockNumber;
    short ID;
    short CodeTable;
    short channels;
    short RepeatNumber;
    short DataLength;
    short FC;
    short Gain;
};

struct DATA_PARSE
{
    unsigned char rec_len;
    unsigned char channle_num;

};
struct DATA_PARSE data_parse;

void ms_delay2(unsigned delay)
{
     unsigned time;
     timer t;
     t :> time;

     for(unsigned int i=0; i<delay; i++)
     {
          time += 100000;
          t when timerafter(time) :> void;
     }
 }
/*******************************************************************************************
 * Function   :unsigned timing(unsigned starttime,unsigned delay)
 * Description:
 * Parameter  ：starttime,delay
 * Return     ：Timer time to return 1,Timer time not to return 0
 * Author     :Zongxiang Xu
 * *****************************************************************************************/
unsigned timing1(unsigned starttime,unsigned us_delay)
{
     unsigned nowtime;
     timer t;
     t :> nowtime;
     //printf("start:%d \n",starttime);
     //printf("time:%d \n",time);
     if(((nowtime-starttime)/100)>=us_delay)
         return 1;
     else
         return 0;
}

void wifi_init()
{
    set_gpio_wifi(WIFI_EN, 0);
    set_gpio_wifi(WIFI_RST, 0);
    ms_delay2(2000);
    set_gpio_wifi(WIFI_EN, 1);
    ms_delay2(500);
    set_gpio_wifi(WIFI_RST, 1);
}

void Uart_Send_Byte(buffered out port:1 wifi_tx, unsigned char byte)
{
    unsigned time;
    timer t;

    t :> time;
    wifi_tx <: 0; /* 开始位 */
    time += BIT_TIME ;
    t when timerafter ( time ) :> void ;

    for ( int i=0; i < DATA_BITS; i ++)  /* 数据位，由低位到高位 */
    {
        wifi_tx <: >> byte;
        time += BIT_TIME ;
        t when timerafter ( time ) :> void ;
    }

    wifi_tx <: 1; /* 停止位 */
    time += BIT_TIME ;
    t when timerafter ( time ) :> void ;
}

unsigned char Uart_Receive_Byte (buffered in port:1 wifi_rx)
{
    timer t;
    unsigned time;
    unsigned char byte = 0;

    wifi_rx when pinseq (0) :> void ; /* 等待开始位 */

    t :> time ;
    time += BIT_TIME /2; /* 增加半个位的延时时间，这样就能与发送端的时序错  开半个位，能更准确地采样数据 */

    for ( int i=0; i <DATA_BITS; i ++) { /* 接收数据 */
        time += BIT_TIME ;
        t when timerafter ( time ) :> void ;
        wifi_rx :> >> byte;
    }

    time += BIT_TIME ; /* 结束位 */
    t when timerafter ( time ) :> void ;
    return byte;
}
void InitUart(buffered out port:1 wifi_tx, buffered in port:1 wifi_rx)
{
    wifi_tx <: 1;
    set_port_pull_none(wifi_rx);
}



void Uart_Send_Frame(buffered out port:1 wifi_tx,unsigned char senddata[],unsigned char len)
{
    for(unsigned i=0;i<len;i++)
    {
        Uart_Send_Byte(wifi_tx, senddata[i]);
    }
}


unsigned char Uart_Receive_Frame(buffered in port:1 wifi_rx,unsigned char receivedata[],unsigned char len)
{
    int i=0;
    unsigned char byte;
    while(len--)
    {
        byte=Uart_Receive_Byte (wifi_rx);
        if(byte!=0)
        {
            receivedata[i]=byte;
            i++;
        }else
        {
            break;
        }
    }
    return i;
}

unsigned char Uart_Protocol_Frame(buffered in port:1 wifi_rx,unsigned char receivedata[])
{
    int i=2;
    unsigned char byte;
    unsigned char len=0;

    receivedata[0]=Uart_Receive_Byte (wifi_rx);
    if(receivedata[0]==0x68)
    {
        receivedata[1]=Uart_Receive_Byte (wifi_rx);
        len=receivedata[1];

        while(len--)
        {
            byte=Uart_Receive_Byte (wifi_rx);
            receivedata[i]=byte;
            i++;
        }
    }
    return i;
}




void Uart_RecData_Frame(buffered in port:1 wifi_rx,unsigned char rec_data[])
{
    int i=0;
    unsigned char byte;
    unsigned char len=10;//数据头长度
    unsigned char datalen=0;


    while(len--)
    {
      byte=Uart_Receive_Byte(wifi_rx);
      if(byte!=0)
      {
          rec_data[i]=byte;
          i++;
      }
    }
    data_parse.channle_num=rec_data[7];

    if(rec_data[2] == 0x2B)
    {
        byte=Uart_Receive_Byte(wifi_rx);
        rec_data[10]=byte;
        if(rec_data[10]!=58)
        {
            byte=Uart_Receive_Byte(wifi_rx);
            rec_data[11]=byte;
            if(rec_data[11]!=58)
            {
                datalen=(rec_data[9]-48)*100+(rec_data[10]-48)*10+(rec_data[11]-47);
                i=i+2;
            }else
            {
                datalen=(rec_data[9]-48)*10+(rec_data[10]-48);
                i=i+2;
            }
        }else
        {
            datalen=rec_data[9]-48;
            i++;
        }

        i=0;
        while(datalen--)
        {
            byte=Uart_Receive_Byte(wifi_rx);
            rec_data[i]=byte;
            i++;
        }
    }
    data_parse.rec_len=i;
}



void Send_Instructions(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx,unsigned char send_instr[],unsigned char rec_len,unsigned char rec_instr[],unsigned int delay_time)
{

    unsigned char send_len;
    unsigned char rec_data[256]={0};

   // while(1)
    {
        send_len=strlen(send_instr);
        printf("\nTX:");
        for(unsigned int i=0;i<send_len;i++)
        {
            printf("%c",send_instr[i]);
        }
        Uart_Send_Frame(wifi_tx,send_instr,send_len);
        rec_len =  Uart_Receive_Frame(wifi_rx,rec_data,rec_len);//Uart_Receive_Frame(wifi_rx);
        ms_delay2(delay_time);
        printf("\nRX:");
        for(unsigned int i=0;i<rec_len;i++)
        {
            printf("%c",rec_data[i]);
        }
       // if(strstr(rec_data,rec_instr))
        {
            //break;
        }
    }
}


void esp8266_ap_client(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx,unsigned char wifi_id_change[])
{
    unsigned char senddata[256]={0};
    unsigned char receivedata[256]={0};
    unsigned char DataLen;


    Send_Instructions(wifi_tx,wifi_rx,"AT\r\n",10,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CWMODE=2\r\n",19,"OK",500);
    Send_Instructions(wifi_tx,wifi_rx,"AT+RST\r\n",14,"OK",2000);

    strncpy(senddata,"AT+CWSAP=\"",10);
    if(wifi_id_change[0]==0x68 && wifi_id_change[2]==0x03)
    {

        strncpy(&senddata[10],&wifi_id_change[5],13);
    }else
    {
        strncpy(&senddata[10],wifi_id,13);
    }

    strncpy(&senddata[23],"\",\"123456789\",11,3\r\n",20);

    DataLen=43;
    printf("\nTX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",senddata[i]);
    }
    Uart_Send_Frame(wifi_tx,senddata,DataLen);

    DataLen=50;
    DataLen =  Uart_Receive_Frame(wifi_rx,receivedata,DataLen);//Uart_Receive_Frame(wifi_rx);
    printf("\nRX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",receivedata[i]);
    }
    ms_delay2(1000);

    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPMUX=0\r\n",19,"OK",1000);

    while(1)
    {
        strncpy(senddata,"AT+CIPSTART=\"TCP\",\"192.168.4.2\",8080\r\n",38);
        DataLen=38;
        printf("\nTX:");
        for(unsigned int i=0;i<DataLen;i++)
        {
            printf("%c",senddata[i]);
        }
        Uart_Send_Frame(wifi_tx,senddata,DataLen);

        DataLen=54;
        DataLen =  Uart_Receive_Frame(wifi_rx,receivedata,DataLen);//Uart_Receive_Frame(wifi_rx);
        printf("\nRX:");
        for(unsigned int i=0;i<DataLen;i++)
        {
            printf("%c",receivedata[i]);
        }
        if(strstr(receivedata,"CONNECT"))
        {
            printf("CONNECT 192.168.4.2:8080 SUCCESS!!\n");
            break;
        }
        ms_delay2(1000);
    }

    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPMODE=1\r\n",20,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSEND\r\n",20,"OK",1000);
}


void esp8266_cipclose(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx)
{
    unsigned char senddata[3]={0};
    unsigned char DataLen;

    for(unsigned char j=0;j<3;j++)
    {
        strncpy(senddata,"+++",3);
        DataLen=3;
        Uart_Send_Frame(wifi_tx,senddata,DataLen);
        printf("\nTX:");
        for(unsigned int i=0;i<DataLen;i++)
        {
            printf("%c",senddata[i]);
        }
        ms_delay2(1000);
    }

    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPCLOSE\r\n",22,"CLOSED",1000);
}
void esp8266_sta_client(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx,unsigned char wifi_id_change[])
{
    unsigned char senddata[256]={0};
    unsigned char receivedata[256]={0};
    unsigned char DataLen;

    Send_Instructions(wifi_tx,wifi_rx,"AT\r\n",10,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CWMODE=3\r\n",19,"OK",500);
    Send_Instructions(wifi_tx,wifi_rx,"AT+RST\r\n",14,"OK",2000);

    strncpy(senddata,"AT+CWSAP=\"",10);
    if(wifi_id_change[0]==0x68 && wifi_id_change[2]==0x03)
    {
        strncpy(&senddata[10],&wifi_id_change[5],13);
    }else
    {
        strncpy(&senddata[10],wifi_id,13);
    }

    strncpy(&senddata[23],"\",\"123456789\",11,3\r\n",20);
    DataLen=43;
    printf("\nTX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",senddata[i]);
    }
    Uart_Send_Frame(wifi_tx,senddata,DataLen);

    DataLen=50;
    DataLen =  Uart_Receive_Frame(wifi_rx,receivedata,DataLen);
    ms_delay2(500);
    printf("\nRX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",receivedata[i]);
    }

    Send_Instructions(wifi_tx,wifi_rx,"AT+CWJAP=\"TP-LINK_3D05\",\"adminadmin\"\r\n",68,"WIFI CONNECT",3000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPMUX=0\r\n",19,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSTART=\"TCP\",\"192.168.1.187\",81\r\n",50,"CONNECT",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPMODE=1\r\n",19,"OK",500);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSEND\r\n",22,"OK",500);
}

void esp8266_sta_ap(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx,unsigned char wifi_id_change[])
{
    unsigned char senddata[256]={0};
    unsigned char receivedata[256]={0};
    unsigned char DataLen;

    Send_Instructions(wifi_tx,wifi_rx,"AT\r\n",10,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CWMODE=3\r\n",19,"OK",500);
    Send_Instructions(wifi_tx,wifi_rx,"AT+RST\r\n",14,"OK",2000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPMUX=1\r\n",19,"OK",1000);

    strncpy(senddata,"AT+CWSAP=\"",10);
    if(wifi_id_change[0]==0x68 && wifi_id_change[2]==0x03)
    {
        strncpy(&senddata[10],&wifi_id_change[5],13);
    }else
    {
        strncpy(&senddata[10],wifi_id,13);
    }

    strncpy(&senddata[23],"\",\"123456789\",11,3\r\n",20);

    DataLen=43;
    printf("\nTX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",senddata[i]);
    }
    Uart_Send_Frame(wifi_tx,senddata,DataLen);


    DataLen=50;
    DataLen =  Uart_Receive_Frame(wifi_rx,receivedata,DataLen);//Uart_Receive_Frame(wifi_rx);
    printf("\nRX:");
    for(unsigned int i=0;i<DataLen;i++)
    {
        printf("%c",receivedata[i]);
    }
    ms_delay2(1000);

    while(1)
    {
        strncpy(senddata,"AT+CWJAP=\"TP-LINK_3D05\",\"adminadmin\"\r\n",38);
        /*strncpy(senddata,"AT+CWSAP=",9);
        strcat(senddata,name);
        strcat(senddata,",");
        strcat(senddata,password);
        strcat(senddata,",11,3\r\n");*/
        DataLen=38;
        printf("\nTX:");
        for(unsigned int i=0;i<DataLen;i++)
        {
            printf("%c",senddata[i]);
        }
        Uart_Send_Frame(wifi_tx,senddata,DataLen);

        DataLen=68;
        DataLen =  Uart_Receive_Frame(wifi_rx,receivedata,DataLen);//Uart_Receive_Frame(wifi_rx);
        ms_delay2(2000);
        printf("\nRX:");
        for(unsigned int i=0;i<DataLen;i++)
        {
            printf("%c",receivedata[i]);
        }
        if(strstr(receivedata,"WIFI CONNECT"))
        {
            break;
        }
    }
    ms_delay2(1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSERVER=1,8080\r\n",26,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSTART=0,\"TCP\",\"192.168.1.187\",81\r\n",55,"OK",1000);
    Send_Instructions(wifi_tx,wifi_rx,"AT+CIPSTO=600\r\n",21,"OK",1000);
}
void Channle_Send(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx,unsigned char channle_num[],unsigned char send_len[])
{
    unsigned char send_instr[20]={0};
    unsigned char rec_instr[30]={0};

    strncpy(send_instr,"AT+CIPSEND=",11);
    strncpy(&send_instr[11],channle_num,1);
    strncpy(&send_instr[12],",",1);
    strcat(&send_instr[13],send_len);

    strncpy(&send_instr[13+strlen(send_len)],"\r\n",2);
   /* printf("\r\n");
    for(unsigned int i=0;i<strlen(send_instr);i++)
    {
        printf("%c",send_instr[i]);
    }*/

    Uart_Send_Frame(wifi_tx,send_instr,strlen(send_instr));

    Uart_Receive_Frame(wifi_rx,rec_instr,(strlen(send_instr)+5));
    /*for(unsigned int i=0;i<strlen(send_instr)+5;i++)
    {
        printf("%c",rec_instr[i]);
    }*/
}


void Double_Protocal_Process(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx)
{
    unsigned char ID_data[256]={0};
    unsigned char send_data[256]={0};
    unsigned char send_len;
    unsigned char rec_data[256]={0};
    unsigned char rec_len;
    static int erase_status=1;
    static int updata_page_count=0;
    struct code_version version;
    unsigned char heart_package[6];
    //memset(receivedata,0,256);
    //read_to_flash(1,receivedata);

    read_flash(16,ID_data);//read ID
    esp8266_sta_ap(wifi_tx,wifi_rx,ID_data);
    set_gpio_led(LED2, 0);

    heart_package[0]=STR_HEAD;
    heart_package[1]=0x04;
    heart_package[2]=HERAT_MSG;
    heart_package[3]=START_TRANS_ACK;
    heart_package[4]=0x00;
    heart_package[5]=0x00;


    version.BlockNumber=0;
    version.ID=0;
    version.CodeTable=113;
    version.channels=5;
    version.RepeatNumber=8;
    version.DataLength=2048;
    version.FC=18000;
    version.Gain=16;

    //esp8266_ap_client(wifi_tx, wifi_rx,ID_data);
    //84 65 31 39 31 31 30 30 30 3130 30 65
    //00 00 00 00 00 71 00 05 00 08 08 00 46 50 00 10
    send_data[0]=STR_HEAD;
    send_data[1]=0x20;
    send_data[2]=REPORT_MSG;
    send_data[3]=START_TRANS;
    send_data[4]=0x00;
    if((ID_data[0]==0x68) && (ID_data[2]==0x03))
    {
        memcpy(&send_data[5],&ID_data[5],13);//SN信息
    }else
    {
        memcpy(&send_data[5],&wifi_id[0],13);//SN信息
    }
    memset(&send_data[18],0x30,5);
    send_data[23]=(unsigned char)(version.CodeTable & 0xff);
    send_data[24]=0x30;
    send_data[25]=(unsigned char)(version.channels & 0xff)+0x30;
    send_data[26]=0x30;
    send_data[27]=(unsigned char)(version.RepeatNumber & 0xff)+0x30;
    send_data[28]=(unsigned char)(version.DataLength >> 8)+0x30;
    send_data[29]=(unsigned char)(version.DataLength & 0xff)+0x30;
    send_data[30]=(unsigned char)(version.FC >> 8);
    send_data[31]=(unsigned char)(version.FC & 0xff);
    send_data[32]=0x30;
    send_data[33]=(unsigned char)(version.Gain & 0xff);

    send_len=34;

    Channle_Send(wifi_tx,wifi_rx,"0","34");



    ms_delay2(1000);
    Uart_Send_Frame(wifi_tx,send_data,send_len);
    printf("\nReport XMOS Message:");
    for(unsigned int i=0;i<send_len;i++)
    {
        printf("%02X ",send_data[i]);
    }
    printf("\n");
    //Uart_Receive_Frame(wifi_rx,send_data,send_len);

    while(1)
    {
        memset(rec_data,0,256);
        //Uart_Receive_Frame(wifi_rx,rec_data,18);
        Uart_RecData_Frame(wifi_rx,rec_data);

        printf("\nchannel:%c len:%d",data_parse.channle_num,data_parse.rec_len);
        printf("\nRX:");
        for(unsigned int i=0;i<data_parse.rec_len;i++)
        {
            printf("%02X ",rec_data[i]);
        }
        printf("\n");

        if(rec_data[0]==0x68)
        {
            switch(rec_data[2])
            {

                case REPORT_MSG://上报设备信息
                    break;
                case HERAT_MSG://心跳

                    if(data_parse.channle_num==48)
                    {
                        Channle_Send(wifi_tx,wifi_rx,"0","6");
                    }else
                    {
                        Channle_Send(wifi_tx,wifi_rx,"1","6");
                    }
                    //Channle_Send(wifi_tx,wifi_rx,"0","6");
                    Uart_Send_Frame(wifi_tx,heart_package,6);
                    printf("\nTX:");
                    for(unsigned int i=0;i<6;i++)
                    {
                        printf("%02X ",heart_package[i]);
                    }
                    printf("\n");
                    break;
                case REBOOT_MSG://复位

                     printf("device_reboot1\n");
                     rec_data[3]=ACTIVE_ACK;
                     if(data_parse.channle_num==48)
                     {
                         Channle_Send(wifi_tx,wifi_rx,"0","6");
                     }else
                     {
                         Channle_Send(wifi_tx,wifi_rx,"1","6");
                     }
                     Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<data_parse.rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                    // esp8266_cipclose(wifi_tx,wifi_rx);
                     ms_delay2(100);
                     device_reboot1();
                     break;
                case GAIN_MODIFY://修改增益,在扇区0的第一页

                     flash_erase_sector(0);
                     write_data_flash(0,rec_data);
                     ms_delay2(100);
                     rec_data[3]=ACTIVE_ACK;
                     if(data_parse.channle_num==48)
                     {
                         Channle_Send(wifi_tx,wifi_rx,"0","6");
                     }else
                     {
                         Channle_Send(wifi_tx,wifi_rx,"1","6");
                     }
                     Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<data_parse.rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                     esp8266_cipclose(wifi_tx,wifi_rx);
                     ms_delay2(100);
                     device_reboot1();
                     break;
                case WIFI_SN://ESN修改，在扇区1的第一页即页的数组16
                    flash_erase_sector(1);
                    write_data_flash(16,rec_data);

                    rec_data[3]=ACTIVE_ACK;
                    if(data_parse.channle_num==48)
                    {
                        Channle_Send(wifi_tx,wifi_rx,"0","18");
                    }else
                    {
                        Channle_Send(wifi_tx,wifi_rx,"1","18");
                    }
                    Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                    printf("\nTX:");
                    for(unsigned int i=0;i<data_parse.rec_len;i++)
                    {
                        printf("%02X ",rec_data[i]);
                    }
                    printf("\n");
                    //esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
                case RESTORE_SET://恢复出厂设置
                    rec_data[3]=ACTIVE_ACK;
                    if(data_parse.channle_num==48)
                    {
                        Channle_Send(wifi_tx,wifi_rx,"0","6");
                    }else
                    {
                        Channle_Send(wifi_tx,wifi_rx,"1","6");
                    }
                    Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                    flash_erase_block(0);
                    ms_delay2(100);
                    flash_erase_block(1);
                    ms_delay2(100);
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    device_reboot1();
                    break;
                case CODE_UPDATA://升级
                    if(rec_data[3]==0x90)
                    {
                        rec_data[3]=0x91;
                        if(data_parse.channle_num==48)
                        {
                            Channle_Send(wifi_tx,wifi_rx,"0","24");
                        }else
                        {
                            Channle_Send(wifi_tx,wifi_rx,"1","24");
                        }
                        Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                        printf("\nTX:");
                        for(unsigned int i=0;i<data_parse.rec_len;i++)
                        {
                            printf("%02X ",rec_data[i]);
                        }
                        printf("\n");
                        flash_erase_block(1);
                        ms_delay2(100);
                        write_data_flash(1*16*16+updata_page_count,rec_data);
                    }else if(rec_data[3]==0x92)
                    {
                        set_gpio_led(LED2, 1);
                       /* if(erase_status)
                        {
                            flash_erase_block(1);
                            ms_delay1(100);
                            write_data_flash(1*16*16+updata_page_count,rec_data);

                            erase_status=0;
                        }else*/
                        {
                            write_data_flash(1*16*16+updata_page_count,rec_data);
                            set_gpio_led(LED2, 0);
                        }
                        rec_data[1]=0x07;
                        rec_data[3]=0x93;
                        rec_data[8]=0x01;
                        rec_len=9;

                        if(data_parse.channle_num==48)
                        {
                            Channle_Send(wifi_tx,wifi_rx,"0","9");
                        }else
                        {
                            Channle_Send(wifi_tx,wifi_rx,"1","9");
                        }

                        Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                        printf("\nTX:");
                        for(unsigned int i=0;i<rec_len;i++)
                        {
                            printf("%02X ",rec_data[i]);
                        }
                        printf("\n");
                        updata_page_count++;


                        }else if(rec_data[3]==0x94)
                        {
                            updata_page_count=0;
                            rec_data[1]=0x06;
                            rec_data[3]=0x95;

                            rec_len=8;

                            if(data_parse.channle_num==48)
                            {
                                Channle_Send(wifi_tx,wifi_rx,"0","8");
                            }else
                            {
                                Channle_Send(wifi_tx,wifi_rx,"1","8");
                            }
                            Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                            printf("\nTX:");
                            for(unsigned int i=0;i<rec_len;i++)
                            {
                                printf("%02X ",rec_data[i]);
                            }
                            printf("\n");
                            esp8266_cipclose(wifi_tx,wifi_rx);
                            ms_delay2(1000);
                            device_reboot1();
                    }

                    break;
                case AMP_SWITCH://功放开关
                     flash_erase_sector(2);
                     write_data_flash(2*16,rec_data);
                     ms_delay2(100);
                     read_flash(2*16,readdata);

                     rec_data[3]=ACTIVE_ACK;
                     //receivedata[5]=0x01;
                     if(data_parse.channle_num==48)
                     {
                         Channle_Send(wifi_tx,wifi_rx,"0","6");
                     }else
                     {
                         Channle_Send(wifi_tx,wifi_rx,"1","6");
                     }
                     Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<data_parse.rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");

                     ms_delay2(100);
                     esp8266_cipclose(wifi_tx,wifi_rx);
                     device_reboot1();
                    break;
                case CONFIG_INQURE://查询参数
                    if(data_parse.channle_num==48)
                    {
                        Channle_Send(wifi_tx,wifi_rx,"0","21");
                    }else
                    {
                        Channle_Send(wifi_tx,wifi_rx,"1","21");
                    }
                     read_flash(1*16*16+0,rec_data);
                     if(rec_data[0]==0xFF)
                     {
                         memset(&rec_data[5],0,5);
                         rec_data[10]=(unsigned char)(version.CodeTable & 0xff);
                         rec_data[11]=0x00;
                         rec_data[12]=(unsigned char)(version.channels & 0xff);
                         rec_data[13]=0x00;
                         rec_data[14]=(unsigned char)(version.RepeatNumber & 0xff);
                         rec_data[15]=(unsigned char)(version.DataLength >> 8);
                         rec_data[16]=(unsigned char)(version.DataLength & 0xff);
                         rec_data[17]=(unsigned char)(version.FC >> 8);
                         rec_data[18]=(unsigned char)(version.FC & 0xff);
                         rec_data[19]=0x00;
                         rec_data[20]=(unsigned char)(version.Gain & 0xff);
                     }
                     rec_data[0]=STR_HEAD;
                     rec_data[1]=0x13;
                     rec_data[2]=CONFIG_INQURE;
                     rec_data[3]=ACTIVE_ACK;
                     rec_data[4]=0x00;
                     Uart_Send_Frame(wifi_tx,rec_data,21);
                     printf("\nTX:");
                     for(unsigned int i=0;i<21;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                    break;
                case SEND_MODE://发送模式
                    rec_data[3]=ACTIVE_ACK;
                    Uart_Send_Frame(wifi_tx,rec_data,data_parse.rec_len);
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
                case 0xC2://重连平台
                    /*rec_data[3]=0x91;
                    Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                    printf("\nTX:");
                    for(unsigned int i=0;i<rec_len;i++)
                    {
                        printf("%02X ",rec_data[i]);
                    }
                    printf("\n");*/
                    //esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
              default:
                  break;
            }
        }
    }
}

void Protocal_Process(buffered out port:1 wifi_tx,buffered in port:1 wifi_rx)
{
    unsigned char ID_data[256]={0};
    unsigned char send_data[256]={0};
    unsigned char send_len;
    unsigned char rec_data[256]={0};
    unsigned char rec_len;
    static int erase_status=1;
    static int updata_page_count=0;
    struct code_version version;
    unsigned char heart_package[6];
    //memset(receivedata,0,256);
    //read_to_flash(1,receivedata);

    read_flash(16,ID_data);//read ID
    esp8266_sta_client(wifi_tx,wifi_rx,ID_data);
    set_gpio_led(LED2, 0);

    heart_package[0]=STR_HEAD;
    heart_package[1]=0x04;
    heart_package[2]=HERAT_MSG;
    heart_package[3]=START_TRANS_ACK;
    heart_package[4]=0x00;
    heart_package[5]=0x00;


    version.BlockNumber=0;
    version.ID=0;
    version.CodeTable=113;
    version.channels=5;
    version.RepeatNumber=8;
    version.DataLength=2048;
    version.FC=18000;
    version.Gain=16;

    //esp8266_ap_client(wifi_tx, wifi_rx,ID_data);
    //84 65 31 39 31 31 30 30 30 3130 30 65
    //00 00 00 00 00 71 00 05 00 08 08 00 46 50 00 10
    send_data[0]=STR_HEAD;
    send_data[1]=0x20;
    send_data[2]=REPORT_MSG;
    send_data[3]=START_TRANS;
    send_data[4]=0x00;
    if((ID_data[0]==0x68) && (ID_data[2]==0x03))
    {
        memcpy(&send_data[5],&ID_data[5],13);//SN信息
    }else
    {
        memcpy(&send_data[5],&wifi_id[0],13);//SN信息
    }
    memset(&send_data[18],0x30,5);
    send_data[23]=(unsigned char)(version.CodeTable & 0xff);
    send_data[24]=0x30;
    send_data[25]=(unsigned char)(version.channels & 0xff)+0x30;
    send_data[26]=0x30;
    send_data[27]=(unsigned char)(version.RepeatNumber & 0xff)+0x30;
    send_data[28]=(unsigned char)(version.DataLength >> 8)+0x30;
    send_data[29]=(unsigned char)(version.DataLength & 0xff)+0x30;
    send_data[30]=(unsigned char)(version.FC >> 8);
    send_data[31]=(unsigned char)(version.FC & 0xff);
    send_data[32]=0x30;
    send_data[33]=(unsigned char)(version.Gain & 0xff);

    send_len=34;


    printf("\nReport XMOS Message:");
    for(unsigned int i=0;i<send_len;i++)
    {
        printf("%02X ",send_data[i]);
    }
    printf("\n");
    Uart_Send_Frame(wifi_tx,send_data,send_len);
    //Uart_Receive_Frame(wifi_rx,send_data,send_len);

    while(1)
    {
        memset(rec_data,0,256);
        rec_len =  Uart_Protocol_Frame(wifi_rx,rec_data);

        printf("\nRX:");
        for(unsigned int i=0;i<rec_len;i++)
        {
            printf("%02X ",rec_data[i]);
        }
        printf("\n");

        if(rec_data[0]==0x68)
        {
            switch(rec_data[2])
            {

                case REPORT_MSG://上报设备信息
                    break;
                case HERAT_MSG://心跳
                    Uart_Send_Frame(wifi_tx,heart_package,6);
                    printf("\nTX:");
                    for(unsigned int i=0;i<6;i++)
                    {
                        printf("%02X ",heart_package[i]);
                    }
                    printf("\n");
                    break;
                case REBOOT_MSG://复位

                     printf("device_reboot1\n");
                     rec_data[3]=ACTIVE_ACK;
                     Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                     esp8266_cipclose(wifi_tx,wifi_rx);
                     ms_delay2(100);
                     device_reboot1();
                     break;
                case GAIN_MODIFY://修改增益,在扇区0的第一页

                     flash_erase_sector(0);
                     write_data_flash(0,rec_data);
                     ms_delay2(100);
                     rec_data[3]=ACTIVE_ACK;
                     Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                     esp8266_cipclose(wifi_tx,wifi_rx);
                     ms_delay2(100);
                     device_reboot1();
                     break;
                case WIFI_SN://ESN修改，在扇区1的第一页即页的数组16
                    flash_erase_sector(1);
                    write_data_flash(16,rec_data);

                    rec_data[3]=ACTIVE_ACK;
                    Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                    printf("\nTX:");
                    for(unsigned int i=0;i<rec_len;i++)
                    {
                        printf("%02X ",rec_data[i]);
                    }
                    printf("\n");
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
                case RESTORE_SET://恢复出厂设置
                    rec_data[3]=ACTIVE_ACK;
                    Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                    flash_erase_block(0);
                    ms_delay2(100);
                    flash_erase_block(1);
                    ms_delay2(100);
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    device_reboot1();
                    break;
                case CODE_UPDATA://升级
                    if(rec_data[3]==0x90)
                    {
                        rec_data[3]=0x91;
                        Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                        printf("\nTX:");
                        for(unsigned int i=0;i<rec_len;i++)
                        {
                            printf("%02X ",rec_data[i]);
                        }
                        printf("\n");
                        flash_erase_block(1);
                        ms_delay2(100);
                        write_data_flash(1*16*16+updata_page_count,rec_data);
                    }else if(rec_data[3]==0x92)
                    {
                        set_gpio_led(LED2, 1);
                       /* if(erase_status)
                        {
                            flash_erase_block(1);
                            ms_delay1(100);
                            write_data_flash(1*16*16+updata_page_count,rec_data);

                            erase_status=0;
                        }else*/
                        {
                            write_data_flash(1*16*16+updata_page_count,rec_data);
                            set_gpio_led(LED2, 0);
                        }
                        rec_data[1]=0x07;
                        rec_data[3]=0x93;
                        rec_data[8]=0x01;
                        rec_len=9;
                        Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                        printf("\nTX:");
                        for(unsigned int i=0;i<rec_len;i++)
                        {
                            printf("%02X ",rec_data[i]);
                        }
                        printf("\n");
                        updata_page_count++;


                        }else if(rec_data[3]==0x94)
                        {
                            updata_page_count=0;
                            rec_data[1]=0x06;
                            rec_data[3]=0x95;

                            rec_len=8;
                            Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                            printf("\nTX:");
                            for(unsigned int i=0;i<rec_len;i++)
                            {
                                printf("%02X ",rec_data[i]);
                            }
                            printf("\n");
                            esp8266_cipclose(wifi_tx,wifi_rx);
                            ms_delay2(1000);
                            device_reboot1();
                    }

                    break;
                case AMP_SWITCH://功放开关
                     flash_erase_sector(2);
                     write_data_flash(2*16,rec_data);
                     ms_delay2(100);
                     read_flash(2*16,readdata);

                     rec_data[3]=ACTIVE_ACK;
                     //receivedata[5]=0x01;
                     Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                     printf("\nTX:");
                     for(unsigned int i=0;i<rec_len;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");

                     ms_delay2(100);
                     esp8266_cipclose(wifi_tx,wifi_rx);
                     device_reboot1();
                    break;
                case CONFIG_INQURE://查询参数

                     read_flash(1*16*16+0,rec_data);
                     if(rec_data[0]==0xFF)
                     {
                         memset(&rec_data[5],0,5);
                         rec_data[10]=(unsigned char)(version.CodeTable & 0xff);
                         rec_data[11]=0x00;
                         rec_data[12]=(unsigned char)(version.channels & 0xff);
                         rec_data[13]=0x00;
                         rec_data[14]=(unsigned char)(version.RepeatNumber & 0xff);
                         rec_data[15]=(unsigned char)(version.DataLength >> 8);
                         rec_data[16]=(unsigned char)(version.DataLength & 0xff);
                         rec_data[17]=(unsigned char)(version.FC >> 8);
                         rec_data[18]=(unsigned char)(version.FC & 0xff);
                         rec_data[19]=0x00;
                         rec_data[20]=(unsigned char)(version.Gain & 0xff);
                     }
                     rec_data[0]=STR_HEAD;
                     rec_data[1]=0x13;
                     rec_data[2]=CONFIG_INQURE;
                     rec_data[3]=ACTIVE_ACK;
                     rec_data[4]=0x00;
                     Uart_Send_Frame(wifi_tx,rec_data,21);
                     printf("\nTX:");
                     for(unsigned int i=0;i<21;i++)
                     {
                         printf("%02X ",rec_data[i]);
                     }
                     printf("\n");
                    break;
                case SEND_MODE://发送模式
                    rec_data[3]=ACTIVE_ACK;
                    Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
                case 0xC2://重连平台
                    /*rec_data[3]=0x91;
                    Uart_Send_Frame(wifi_tx,rec_data,rec_len);
                    printf("\nTX:");
                    for(unsigned int i=0;i<rec_len;i++)
                    {
                        printf("%02X ",rec_data[i]);
                    }
                    printf("\n");*/
                    esp8266_cipclose(wifi_tx,wifi_rx);
                    ms_delay2(100);
                    device_reboot1();
                    break;
              default:
                  break;
            }
        }
    }
}


