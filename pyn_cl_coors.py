import io_formats as io

class cl_traj():
    def __init__(self,file_input=None,frame=None,begin=None,end=None,increment=1,units=None,verbose=True):

        self.io_file=None
        self.io_opened=False
        self.io_end=False
        self.io_err=False
        self.io_pos=None
        self.io_vars=[]

        self.name=None
        self.type=None
        self.title=None
        self.num_frames=0
        self.precision=None
        self.frame=[]

        if file_input!=None:
            self.name=file_input
            self.type=file_input.split('.')[1]
            self.open()

        if frame!=None or begin!=None or end!=None:
            self.upload_frame(frame,begin,end,increment,units)

        if verbose:
            self.info()
            

    def info(self,index=None):

        if index==None:
            print '#',self.num_frames,'frames/models loaded.'
        else:
            print '#',self.num_frames,'frames/models in traj',index
            

    def open(self):
        
        self.io_file,self.io_vars,self.io_pos,self.io_err=getattr(io,'coor_'+self.type).open_traj(self.name)
        if self.io_err: print '# Error opening the file'; return
        self.io_opened=True


    #def delete_frame(self,frame='ALL',begin=None,end=None,increment=1,units=None):
    # 
    #    if frame in ['all','All','ALL']:
    #        del(self.frame)
    #        self.frame=[]
    #        self.num_frames=0
    #        return
    # 
    #    elif type(frame) in [list,tuple]:
    #        for ii in frame:
    #            self.frame.__delitem__(ii)
    #            self.num_frames-=1
    #        return
    # 
    #    pass

    def reload_frame(self,new='next',old=0):

        if new=='next':
            temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos)
            if self.io_err: 
                print '# Error reading the file'
                return
            if self.io_end: 
                print '# End of file'
                self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                if self.io_err: return '# Error closing file'
                self.io_opened=False
                self.io_file=None
                return
            self.frame[old]=temp_frame
            del(temp_frame)
        else:
            print '# Not supported yet.'
            return

    def upload_frame(self,frame='next',begin=None,end=None,increment=1,units=None):

        if type(frame) in [int]: frame=[frame]

        if begin!=None or end!=None:
            if units==None:
                print "# Option units=None."
                print "# Choose units in ['ns','ps','md_steps','frames']"
                pass
            print "# Not supported yet."
            pass

        else:
            ### Uploading next frame
            if frame in ['next','Next','NEXT']:
                temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_next(self.io_file,self.io_vars,self.io_pos)
                if self.io_err: 
                    print '# Error reading the file'
                    return
                if self.io_end: 
                    print '# End of file'
                    self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                    if self.io_err: return '# Error closing file'
                    self.io_opened=False
                    self.io_file=None
                    return
                self.frame.append(temp_frame)
                self.num_frames+=1
            ### Uploading all frames
            elif frame in ['all','All','ALL']:
                temp_frames,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_all(self.io_file,self.io_vars,self.io_pos)
                if self.io_err: return '# Error reading file'
                self.io_err=getattr(io,'coor_'+self.type).close_traj(self.io_file)
                if self.io_err: return '# Error closing file'
                self.io_opened=False
                self.io_file=None
                self.num_frames+=len(temp_frames)
                for ii in temp_frames:
                    self.frame.append(ii)
                self.precision=self.frame[0].precision
                del(temp_frames)
            ### Uploading a list of frames
            elif type(frame) in [tuple,list]:
                for ii in frame:
                    temp_frame,self.io_pos,self.io_err,self.io_end=getattr(io,'coor_'+self.type).read_frame(self.io_file,ii,self.io_vars,self.io_pos)
                    if self.io_err: return '# Error reading file'
                    self.frame.append(temp_frame)
                    self.num_frames+=1

            else:
                print "# Options not supported yet."
                pass

        pass



        
